#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Compute directional spatial Granger causality on 2D grid data using symbolized entropy measures.
 *
 * This function estimates the directional spatial Granger causality between two variables `x` and `y`,
 * each defined on a 2D grid (i.e., spatial lattice), based on symbolic transfer entropy principles.
 * It computes causality in both directions (X → Y and Y → X) by comparing changes in entropy and
 * conditional entropy across spatially-embedded versions of the input grids.
 *
 * The method includes:
 * - Symbolizing the original and embedded grids using local neighborhoods (controlled by `k`).
 * - Computing joint and marginal entropies involving the original grids and their local embeddings.
 * - Deriving directional causality using entropy difference formulas:
 *   - SC_{x→y} = [H(y, wy) − H(wy)] − [H(y, wy, wx) − H(wy, wx)]
 *   - SC_{y→x} = [H(x, wx) − H(wx)] − [H(x, wx, wy) − H(wx, wy)]
 *     where `wx` and `wy` are spatial embeddings (local lag structures) of `x` and `y`.
 *
 * Symbolization helps convert continuous grid values into discrete symbolic patterns, enabling
 * robust estimation of information-theoretic metrics under spatial autocorrelation.
 *
 * @param x         2D grid of variable X, represented as a vector of rows.
 * @param y         2D grid of variable Y, same size as x.
 * @param k         Neighborhood size for symbolization (e.g., k = 1 implies 3×3 window).
 * @param base      Logarithm base used for entropy calculation (default is 2 for bits).
 * @param symbolize Whether to apply symbolization (currently assumed true by default).
 *
 * @return A std::vector<double> of two values:
 *         - sc_x_to_y: Directional spatial Granger causality from X to Y.
 *         - sc_y_to_x: Directional spatial Granger causality from Y to X.
 */
std::vector<double> SCTSingle4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    size_t k,
    double base = 2,
    bool symbolize = true
) {
  size_t rows = x.size();
  size_t cols = x[0].size();

  std::vector<double> wx(rows*cols);
  std::vector<std::vector<double>> Ex = GenGridEmbeddings(x,1,1);
  for (const auto& row : Ex) {
    wx.insert(wx.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> xw = GridVec2Mat(wx,rows);

  std::vector<double> wy(rows*cols);
  std::vector<std::vector<double>> Ey = GenGridEmbeddings(y,1,1);
  for (const auto& row : Ey) {
    wy.insert(wy.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> yw = GridVec2Mat(wy,rows);

  std::vector<double> sx = GenGridSymbolization(x, k);
  std::vector<double> sy = GenGridSymbolization(y, k);
  std::vector<double> swx = GenGridSymbolization(xw, k);
  std::vector<double> swy = GenGridSymbolization(yw, k);

  std::vector<std::vector<double>> sp_series(rows*cols,std::vector<double>(4));
  for (size_t i = 0; i < x.size(); ++i){
    sp_series[i] = {sx[i],sy[i],swx[i],swy[i]}; // 0:x 1:y 2:wx 3:wy
  }

  double Hwx, Hwy, Hxwx, Hywy, Hwxwy, Hwxwyx, Hwxwyy;
  Hxwx = CppJoinEntropy_Disc(sp_series, {0,2}, base, false); // H(x,wx)
  Hywy = CppJoinEntropy_Disc(sp_series, {1,3}, base, false); // H(y,wy)
  Hwx = CppEntropy_Disc(swx, base, false); // H(wx)
  Hwy = CppEntropy_Disc(swy, base, false); // H(wy)
  Hwxwy = CppJoinEntropy_Disc(sp_series,{2,3}, base, false); // H(wx,wy)
  Hwxwyx = CppJoinEntropy_Disc(sp_series,{0,2,3}, base, false); // H(wx,wy,x)
  Hwxwyy = CppJoinEntropy_Disc(sp_series,{1,2,3}, base, false); // H(wx,wy,y)

  double sc_x_to_y = (Hywy - Hwy) - (Hwxwyy - Hwxwy);
  double sc_y_to_x = (Hxwx - Hwx) - (Hwxwyx - Hwxwy);

  return {sc_x_to_y,sc_y_to_x};
}

/**
 * @brief Compute spatial Granger causality for gridded data using spatial block bootstrap.
 *
 * This function performs statistical inference on spatial Granger causality between two gridded
 * variables `x` and `y`, using a symbolic entropy-based method. It estimates the causality measures
 * from x to y and from y to x, and evaluates their statistical significance via block-based
 * spatial bootstrap.
 *
 * The method includes:
 * - Computing the true causality statistic from the original gridded data.
 * - Performing block-wise spatial bootstrap to generate bootstrapped versions of the data.
 * - Estimating the causality statistic for each bootstrap sample.
 * - Computing the empirical p-values by comparing the bootstrapped statistics with the true statistic.
 *
 * The spatial bootstrap preserves local spatial structure using a provided spatial block assignment vector.
 * Grid values are flattened, shuffled in blocks, and reshaped to their original grid layout for computation.
 *
 * The function returns the estimated spatial Granger causality statistics and their bootstrap p-values:
 * - `sc_x_to_y`: Estimated causality from X to Y.
 * - `p_x_to_y`: Bootstrap p-value for X → Y.
 * - `sc_y_to_x`: Estimated causality from Y to X.
 * - `p_y_to_x`: Bootstrap p-value for Y → X.
 *
 * This implementation allows parallel computation for the bootstrap estimation and can display a progress bar
 * to track the process.
 *
 * @param x           2D grid (matrix) of variable X.
 * @param y           2D grid (matrix) of variable Y, same size as x.
 * @param block       Vector assigning each grid cell to a spatial block for bootstrapping.
 * @param k           Neighborhood window size used for symbolization (typically 3 or 5).
 * @param threads     Number of threads to use for parallel bootstrap estimation.
 * @param boot        Number of bootstrap iterations (default: 399).
 * @param base        Base of the logarithm used in entropy computation (default: 2 for bits).
 * @param seed        Seed for the random number generator to ensure reproducibility (default: 42).
 * @param progressbar Whether to show a progress bar during bootstrap computation (default: true).
 *
 * @return A vector of four values:
 *         - sc_x_to_y: Estimated causality from X to Y.
 *         - p_x_to_y: Bootstrap p-value for X → Y.
 *         - sc_y_to_x: Estimated causality from Y to X.
 *         - p_y_to_x: Bootstrap p-value for Y → X.
 */
std::vector<double> SCT4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    const std::vector<int>& block,
    int k,
    int threads,
    int boot = 399,
    double base = 2,
    unsigned int seed = 42,
    bool progressbar = true
){
  // Initialize the bootstrapped realizations of the spatial granger causality statistic
  std::vector<std::vector<double>> sc_bootstraps(boot);

  const size_t rows = x.size();
  const size_t cols = x[0].size();
  auto monte_boots = [&](int n){
    // Use different seed for each iteration to ensure different random samples
    unsigned int current_seed = seed + n;
    // Generate a spatial block bootstrap resample of indices
    std::vector<int> boot_indice = SpatialBlockBootstrap(block,current_seed);
    // Obtain the bootstrapped realization series
    std::vector<double> x_bs(boot_indice.size());
    std::vector<double> y_bs(boot_indice.size());
    for (size_t i = 0; i < boot_indice.size(); ++i){
      std::vector<int> cellindice = RowColFromGrid(i,cols);
      x_bs[i] = x[cellindice[0]][cellindice[1]];
      y_bs[i] = y[cellindice[0]][cellindice[1]];
    }
    std::vector<std::vector<double>> x_boot = GridVec2Mat(x_bs,static_cast<int>(rows));
    std::vector<std::vector<double>> y_boot = GridVec2Mat(y_bs,static_cast<int>(rows));
    // Estimate the bootstrapped realization of the spatial granger causality statistic
    sc_bootstraps[n] = SCTSingle4Grid(x_boot,y_boot,static_cast<size_t>(std::abs(k)),base);
  };

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Parallel computation with or without a progress bar
  if (progressbar) {
    RcppThread::ProgressBar bar(boot, 1);
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
    }, threads_sizet);
  }

  // The "true" spatial granger causality statistic
  std::vector<double> sc = SCTSingle4Grid(x,y,static_cast<size_t>(std::abs(k)),base);
  double scx = sc[0];
  double scy = sc[1];
  // Compute the estimated bootstrap p–value
  double b_xy = 0;
  double b_yx = 0;
  for (size_t i = 0; i < sc_bootstraps.size(); ++i){
    if (sc_bootstraps[i][0] > scx){
      b_xy += 1;
    }
    if (sc_bootstraps[i][1] > scy) {
      b_yx += 1;
    }
  }

  return {scx,b_xy / boot,scy,b_yx / boot};
}

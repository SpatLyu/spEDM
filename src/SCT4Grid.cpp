#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Compute directional spatial Granger causality for 2D grid data using symbolic entropy measures.
 *
 * This function estimates the bidirectional spatial Granger causality between two spatial variables
 * `x` and `y` observed on a 2D lattice/grid. It does so via a symbolic approximation of transfer entropy,
 * which evaluates whether the spatial neighborhood of one variable improves the prediction of the other.
 *
 * The core procedure includes:
 * - Constructing spatial embeddings (`wx`, `wy`) of `x` and `y` using a local neighborhood structure.
 * - Optionally symbolizing all grids to enhance robustness under spatial autocorrelation.
 * - Computing marginal and joint (conditional) entropies using the symbolized or raw values.
 * - Calculating directional Granger causality based on information gain:
 *
 *   Causality from X to Y:
 *     SC_{x→y} = [H(y, wy) − H(wy)] − [H(y, wy, wx) − H(wy, wx)]
 *
 *   Causality from Y to X:
 *     SC_{y→x} = [H(x, wx) − H(wx)] − [H(x, wx, wy) − H(wx, wy)]
 *
 * where `wx` and `wy` are spatial embeddings (i.e., local neighborhood vectors) of `x` and `y`.
 * These terms quantify how much the inclusion of one variable’s spatial context helps predict the other.
 *
 * @param x         2D grid (matrix) representing variable X.
 * @param y         2D grid (matrix) representing variable Y.
 * @param k         Embedding neighborhood radius (e.g., k = 1 means 3×3 window).
 * @param base      Logarithm base for entropy (default is 2, for units in bits).
 * @param symbolize If true, discretize values using symbolic transformation before entropy computation.
 *
 * @return A std::vector<double> of two values:
 *         - sc_x_to_y: spatial Granger causality from X to Y.
 *         - sc_y_to_x: spatial Granger causality from Y to X.
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

  std::vector<double> wx;
  std::vector<std::vector<double>> Ex = GenGridEmbeddings(x,1,1);
  for (const auto& row : Ex) {
    wx.insert(wx.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> xw = GridVec2Mat(wx,rows);

  std::vector<double> wy;
  std::vector<std::vector<double>> Ey = GenGridEmbeddings(y,1,1);
  for (const auto& row : Ey) {
    wy.insert(wy.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> yw = GridVec2Mat(wy,rows);

  double Hwx, Hwy, Hxwx, Hywy, Hwxwy, Hwxwyx, Hwxwyy;

  if (symbolize){
    std::vector<double> sx = GenGridSymbolization(x, k);
    std::vector<double> sy = GenGridSymbolization(y, k);
    std::vector<double> swx = GenGridSymbolization(xw, k);
    std::vector<double> swy = GenGridSymbolization(yw, k);

    std::vector<std::vector<double>> sp_series(rows*cols,std::vector<double>(4));
    for (size_t i = 0; i < x.size(); ++i){
      sp_series[i] = {sx[i],sy[i],swx[i],swy[i]}; // 0:x 1:y 2:wx 3:wy
    }

    Hxwx = CppJoinEntropy_Disc(sp_series, {0,2}, base, true); // H(x,wx)
    Hywy = CppJoinEntropy_Disc(sp_series, {1,3}, base, true); // H(y,wy)
    Hwx = CppEntropy_Disc(swx, base, true); // H(wx)
    Hwy = CppEntropy_Disc(swy, base, true); // H(wy)
    Hwxwy = CppJoinEntropy_Disc(sp_series,{2,3}, base, true); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Disc(sp_series,{0,2,3}, base, true); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Disc(sp_series,{1,2,3}, base, true); // H(wx,wy,y)
  } else {
    std::vector<std::vector<double>> sp_series(rows*cols,std::vector<double>(4));
    for (size_t i = 0; i < x.size(); ++i){
      int ri = i / cols;
      int ci = i % cols;
      sp_series[i] = {x[ri][ci],y[ri][ci],wx[i],wy[i]}; // 0:x 1:y 2:wx 3:wy
    }

    Hxwx = CppJoinEntropy_Cont(sp_series, {0,2}, k, base, true); // H(x,wx)
    Hywy = CppJoinEntropy_Cont(sp_series, {1,3}, k, base, true); // H(y,wy)
    Hwx = CppEntropy_Cont(wx, k, base, true); // H(wx)
    Hwy = CppEntropy_Cont(wy, k, base, true); // H(wy)
    Hwxwy = CppJoinEntropy_Cont(sp_series,{2,3}, k, base, true); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Cont(sp_series,{0,2,3}, k, base, true); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Cont(sp_series,{1,2,3}, k, base, true); // H(wx,wy,y)
  }

  // transformed to fall within [-1, 1]
  double sc_x_to_y = ((Hywy - Hwy) - (Hwxwyy - Hwxwy)) / (Hywy - Hwy);
  double sc_y_to_x = ((Hxwx - Hwx) - (Hwxwyx - Hwxwy)) / ((Hxwx - Hwx));
  // double sc_x_to_y = (Hywy - Hwy) - (Hwxwyy - Hwxwy);
  // double sc_y_to_x = (Hxwx - Hwx) - (Hwxwyx - Hwxwy);

  return {sc_x_to_y,sc_y_to_x};
}

/**
 * @brief Compute spatial Granger causality for gridded data using spatial block bootstrap.
 *
 * This function estimates the directional spatial Granger causality between two gridded variables `x` and `y`,
 * by applying a symbolic entropy-based method, and assesses the statistical significance of the causality using
 * spatial block bootstrap techniques. It calculates the causality in both directions: X → Y and Y → X.
 * Additionally, the function evaluates the significance of the estimated causality statistics by comparing them
 * to bootstrap realizations of the causality.
 *
 * The method involves the following steps:
 * - **Computation of true causality**: The function first calculates the spatial Granger causality statistic
 *   using the original data grids `x` and `y`.
 * - **Spatial block bootstrap resampling**: The grid values are resampled with spatial block bootstrapping.
 *   Each resample preserves local spatial structure and generates new bootstrap realizations of the causality statistic.
 * - **Estimation of causality for bootstrapped samples**: The causality statistic is estimated for each of the
 *   bootstrapped realizations, which involves calculating the symbolic entropy measures and their differences.
 * - **Empirical p-values**: The final p-values for both directional causality estimates (X → Y and Y → X) are
 *   derived by comparing the bootstrapped statistics with the true causality statistics.
 *
 * This approach takes into account spatial autocorrelation and allows the use of parallel processing for faster
 * bootstrap estimation. The spatial bootstrap method involves reshuffling grid cells into spatial blocks,
 * preserving local dependencies, and calculating causality for each realization.
 *
 * @param x           2D grid (matrix) of variable X.
 * @param y           2D grid (matrix) of variable Y, same size as x.
 * @param block       Vector assigning each grid cell to a spatial block for bootstrapping.
 * @param k           Neighborhood window size used for symbolization (typically 3 or 5).
 * @param threads     Number of threads to use for parallel bootstrap estimation.
 * @param boot        Number of bootstrap iterations (default: 399).
 * @param base        Base of the logarithm used in entropy computation (default: 2 for bits).
 * @param seed        Seed for the random number generator to ensure reproducibility (default: 42).
 * @param symbolize   Whether to use symbolic transformation for the grids (default: true).
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
    bool symbolize = true,
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
    sc_bootstraps[n] = SCTSingle4Grid(x_boot,y_boot,static_cast<size_t>(std::abs(k)),base,symbolize);
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
  std::vector<double> sc = SCTSingle4Grid(x,y,static_cast<size_t>(std::abs(k)),base,symbolize);
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

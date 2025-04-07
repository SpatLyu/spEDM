#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Compute the spatial Granger causality statistic between two 2D grid-structured variables.
 *
 * This function estimates the directional influence between two spatial variables `x` and `y` defined
 * on a 2D grid (matrix form), using symbolization-based entropy and conditional entropy measures.
 * The approach is based on symbolic transfer entropy, where symbolization is applied locally to
 * encode neighborhood patterns, followed by computation of information-theoretic quantities.
 *
 * Specifically, it calculates:
 * - H(X), H(Y): Entropy of the symbolized versions of x and y.
 * - H(Y|X), H(X|Y): Conditional entropy between x and y.
 * - SC_{x→y} = H(Y) - H(Y|X): Directional causality from x to y.
 * - SC_{y→x} = H(X) - H(X|Y): Directional causality from y to x.
 *
 * The function allows for the option of symbolization of the input grids, which transforms the grid values
 * into symbolic representations for better capturing spatial neighborhood patterns.
 * The function returns the directional Granger causality values, one for each direction (from x to y and from y to x).
 *
 * @param x     2D grid of variable X (represented as a vector of rows).
 * @param y     2D grid of variable Y (same size as x).
 * @param k     Neighborhood window size for local symbolization (default is 1).
 * @param base  Logarithm base used in entropy computation (default is 2 for bits).
 * @param symbolize Flag indicating whether symbolization of the input grids is applied (default is true).
 *
 * @return A vector with two values:
 *         - sc_x_to_y: Spatial Granger causality from X to Y.
 *         - sc_y_to_x: Spatial Granger causality from Y to X.
 */
std::vector<double> SCTSingle4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    size_t k,
    double base = 2,
    bool symbolize = true
) {
  std::vector<double> sx, sy;
  if (symbolize) {
    sx = GenGridSymbolization(x, k);
    sy = GenGridSymbolization(y, k);
  } else {
    sx = GridMat2Vec(x);
    sy = GridMat2Vec(y);
  }

  double Hx = CppEntropy_Disc(sx, base, false); // H(x)
  double Hy = CppEntropy_Disc(sy, base, false); // H(y)
  double Hxy = CppConditionalEntropy_Disc(sx, sy, base, false); // H(x | y)
  double Hyx = CppConditionalEntropy_Disc(sy, sx, base, false); // H(y | x)

  double sc_x_to_y = Hy - Hyx;
  double sc_y_to_x = Hx - Hxy;

  return {sc_x_to_y, sc_y_to_x};
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
 * @param symbolize   Flag indicating whether symbolization of the input grids is applied (default: true).
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

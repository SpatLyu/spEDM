#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the directional symbolic causality strength between two spatial variables
 *        based on lattice neighborhood structure and symbolization.
 *
 * This function implements a symbolic causality test statistic for spatial data over lattice structures.
 * It follows a non-parametric entropy-based approach to evaluate whether variable x Granger-causes y,
 * as defined in the referenced paper using symbolized spatial neighbors.
 *
 * Specifically, the function calculates:
 *   - H(y): the entropy of the symbolized y values.
 *   - H(y | x): the conditional entropy of y given x.
 *   - H(x): the entropy of the symbolized x values.
 *   - H(x | y): the conditional entropy of x given y.
 * Then it computes:
 *   - sc_x_to_y = H(y) - H(y | x): symbolic causality strength from x to y.
 *   - sc_y_to_x = H(x) - H(x | y): symbolic causality strength from y to x.
 *
 * These correspond to the delta entropy test statistic described in Equation (24)
 * of the attached reference, where a significant positive value of sc_x_to_y indicates
 * that x contains additional information about y not present in y’s own spatial lags.
 *
 * @param x    Input vector of variable x, assumed to be on a lattice/grid.
 * @param y    Input vector of variable y, assumed to be on the same lattice/grid as x.
 * @param nb   Neighborhood list for each lattice point (e.g., rook or queen adjacency).
 * @param k    The number of neighbors considered in the lattice symbolization.
 * @param base The logarithm base used in entropy calculation (default is 2).
 *
 * @return A std::vector<double> containing two values:
 *         [0] sc_x_to_y — symbolic causality strength from x to y;
 *         [1] sc_y_to_x — symbolic causality strength from y to x.
 *
 */
std::vector<double> SCTSingle4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    double base = 2
){
  std::vector<double> sx = GenLatticeSymbolization(x,nb,k);
  std::vector<double> sy = GenLatticeSymbolization(y,nb,k);
  double Hx = CppEntropy_Disc(sx,base,false); // H(x)
  double Hy = CppEntropy_Disc(sy,base,false); // H(y)
  double Hxy = CppConditionalEntropy_Disc(sx, sy, base, false); // H(x | y)
  double Hyx = CppConditionalEntropy_Disc(sy, sx, base, false); // H(y | x)
  double sc_x_to_y = Hy - Hyx;
  double sc_y_to_x = Hx - Hxy;

  return {sc_x_to_y,sc_y_to_x};
}

/**
 * @brief Perform spatial Granger causality testing on lattice data with bootstrap-based significance evaluation.
 *
 * This function computes the spatial Granger causality statistic between two variables `x` and `y`
 * over a lattice structure defined by neighborhood relationships. To evaluate the statistical
 * significance, it applies a spatial block bootstrap procedure to generate empirical p-values
 * under the null hypothesis of no causality.
 *
 * The bootstrap procedure respects the spatial contiguity structure provided by the `block` vector,
 * where each unique block ID defines a spatially contiguous subregion.
 * For each bootstrap iteration, a new realization is generated using spatial block bootstrap,
 * and the test statistic is recalculated.
 *
 * The final output includes the original causality statistics and their associated bootstrap p-values.
 *
 * @param x           Input variable X (time series or spatial measurements).
 * @param y           Input variable Y (potentially Granger-caused by X).
 * @param nb          Neighborhood structure (adjacency list) of the lattice.
 * @param block       Block IDs for spatial block bootstrap (same length as x/y).
 * @param k           The number of neighbors considered in the lattice symbolization.
 * @param threads     Number of threads to use for parallel computation.
 * @param boot        Number of bootstrap replicates to perform (default is 399).
 * @param base        Base for temporal lag scaling (default is 2).
 * @param seed        Random seed for reproducibility (default is 42).
 * @param progressbar Whether to display a progress bar (default is true).
 *
 * @return A vector of 4 values:
 *         - scx: Observed statistic for X → Y causality.
 *         - p_xy: Empirical p-value for X → Y based on bootstrap.
 *         - scy: Observed statistic for Y → X causality.
 *         - p_yx: Empirical p-value for Y → X based on bootstrap.
 */
std::vector<double> SCT4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& block,
    size_t k,
    size_t threads,
    int boot = 399,
    double base = 2,
    unsigned int seed = 42,
    bool progressbar = true
){
  // Initialize the bootstrapped realizations of the spatial granger causality statistic
  std::vector<std::vector<double>> sc_bootstraps(boot);

  auto monte_boots = [&](int n){
    // Use different seed for each iteration to ensure different random samples
    unsigned int current_seed = seed + n;
    // Generate a spatial block bootstrap resample of indices
    std::vector<int> boot_indice = SpatialBlockBootstrap(block,current_seed);
    // Obtain the bootstrapped realization series
    std::vector<double> x_boot(x.size());
    std::vector<double> y_boot(y.size());
    for (size_t i = 0; i < boot_indice.size(); ++i){
      x_boot[i] = x[boot_indice[i]];
      y_boot[i] = y[boot_indice[i]];
    }
    // Estimate the bootstrapped realization of the spatial granger causality statistic
    sc_bootstraps[n] = SCTSingle4Lattice(x_boot,y_boot,nb,k,base);
  };

  // Parallel computation with or without a progress bar
  if (progressbar) {
    RcppThread::ProgressBar bar(boot, 1);
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
      bar++;
    }, threads);
  } else {
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
    }, threads);
  }

  // The "true" spatial granger causality statistic
  std::vector<double> sc = SCTSingle4Lattice(x,y,nb,k,base);
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

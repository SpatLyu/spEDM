#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the directional symbolic causality strength between two spatial variables
 *        over a lattice structure using lattice-based embedding and optional symbolization.
 *
 * This function implements a symbolic causality analysis between two spatial variables `x` and `y`,
 * aligned on a regular lattice (e.g., grid). It is designed for spatial processes and uses
 * lattice-based neighbor embedding and symbolization to quantify directional influence.
 *
 * The method proceeds as follows:
 * 1. **Embedding**: For each location, the function embeds `x` and `y` into a new series `wx` and `wy`,
 *    which incorporate spatial lagged information using a neighborhood structure (`nb`).
 * 2. **Symbolization (optional)**: If `symbolize = true`, both the original (`x`, `y`) and embedded (`wx`, `wy`)
 *    vectors are symbolized using lattice-based quantization with `k` symbols (bins).
 * 3. **Entropy & Joint Entropy Calculation**:
 *    - \( H(wx), H(wy), H(wx, wy) \)
 *    - \( H(x, wx), H(y, wy) \)
 *    - \( H(wx, wy, x), H(wx, wy, y) \)
 * 4. **Symbolic Causality Strength Computation**:
 *    - \( sc_{x \rightarrow y} = [H(y, wy) - H(wy)] - [H(wx, wy, y) - H(wx, wy)] \)
 *    - \( sc_{y \rightarrow x} = [H(x, wx) - H(wx)] - [H(wx, wy, x) - H(wx, wy)] \)
 *    These represent the reduction in uncertainty (entropy) of `y` (or `x`) due to knowledge of `x` (or `y`)
 *    in the context of spatial lags.
 *
 * @param x         Input vector of spatial variable x, aligned with a lattice/grid.
 * @param y         Input vector of spatial variable y, aligned with the same lattice/grid.
 * @param nb        Neighborhood list for each lattice location (e.g., rook or queen neighbors).
 * @param k         Number of discrete symbols used in lattice symbolization.
 * @param base      Base of logarithm used in entropy computation (default is 2).
 * @param symbolize If true, apply lattice-based symbolization before entropy calculation.
 *
 * @return A std::vector<double> of size 2:
 *         - [0] Symbolic causality strength from x to y (sc_x_to_y)
 *         - [1] Symbolic causality strength from y to x (sc_y_to_x)
 */
std::vector<double> SCTSingle4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    double base = 2,
    bool symbolize = true
){
  std::vector<double> wx(x.size());
  std::vector<std::vector<double>> Ex = GenLatticeEmbeddings(x,nb,1,1);
  for (const auto& row : Ex) {
    wx.insert(wx.end(), row.begin(), row.end());
  }

  std::vector<double> wy(y.size());
  std::vector<std::vector<double>> Ey = GenLatticeEmbeddings(y,nb,1,1);
  for (const auto& row : Ey) {
    wy.insert(wy.end(), row.begin(), row.end());
  }

  double Hwx, Hwy, Hxwx, Hywy, Hwxwy, Hwxwyx, Hwxwyy;
  if (symbolize){
    std::vector<double> sx = GenLatticeSymbolization(x,nb,k);
    std::vector<double> sy = GenLatticeSymbolization(y,nb,k);
    std::vector<double> swx = GenLatticeSymbolization(wx,nb,k);
    std::vector<double> swy = GenLatticeSymbolization(wy,nb,k);

    std::vector<std::vector<double>> sp_series(x.size(),std::vector<double>(4));
    for (size_t i = 0; i < x.size(); ++i){
      sp_series[i] = {sx[i],sy[i],swx[i],swy[i]}; // 0:x 1:y 2:wx 3:wy
    }

    Hxwx = CppJoinEntropy_Disc(sp_series, {0,2}, base, false); // H(x,wx)
    Hywy = CppJoinEntropy_Disc(sp_series, {1,3}, base, false); // H(y,wy)
    Hwx = CppEntropy_Disc(swx, base, false); // H(wx)
    Hwy = CppEntropy_Disc(swy, base, false); // H(wy)
    Hwxwy = CppJoinEntropy_Disc(sp_series,{2,3}, base, false); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Disc(sp_series,{0,2,3}, base, false); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Disc(sp_series,{1,2,3}, base, false); // H(wx,wy,y)
  } else {
    std::vector<std::vector<double>> sp_series(x.size(),std::vector<double>(4));
    for (size_t i = 0; i < x.size(); ++i){
      sp_series[i] = {x[i],y[i],wx[i],wy[i]}; // 0:x 1:y 2:wx 3:wy
    }

    Hxwx = CppJoinEntropy_Disc(sp_series, {0,2}, base, false); // H(x,wx)
    Hywy = CppJoinEntropy_Disc(sp_series, {1,3}, base, false); // H(y,wy)
    Hwx = CppEntropy_Disc(wx, base, false); // H(wx)
    Hwy = CppEntropy_Disc(wy, base, false); // H(wy)
    Hwxwy = CppJoinEntropy_Disc(sp_series,{2,3}, base, false); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Disc(sp_series,{0,2,3}, base, false); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Disc(sp_series,{1,2,3}, base, false); // H(wx,wy,y)
  }
  double sc_x_to_y = (Hywy - Hwy) - (Hwxwyy - Hwxwy);
  double sc_y_to_x = (Hxwx - Hwx) - (Hwxwyx - Hwxwy);

  return {sc_x_to_y,sc_y_to_x};
}

/**
 * @brief Perform spatial Granger causality testing on lattice data with bootstrap-based significance evaluation.
 *
 * This function computes the spatial Granger causality statistic between two variables, `x` and `y`,
 * over a lattice structure defined by neighborhood relationships. The statistical significance of the causality
 * is evaluated through a spatial block bootstrap procedure, generating empirical p-values under the null hypothesis
 * of no causality.
 *
 * The spatial block bootstrap preserves the neighborhood relationships defined by the `block` vector, where each
 * unique block ID corresponds to a spatially contiguous subregion. For each bootstrap iteration, a new realization
 * is generated using this procedure, and the test statistic is recalculated to assess the significance of causality.
 *
 * The final output includes the observed causality statistics for both directions (X → Y and Y → X) and their associated
 * bootstrap p-values, which indicate the likelihood of observing such causality by random chance.
 *
 * @param x           Input variable X (time series or spatial measurements).
 * @param y           Input variable Y (potentially Granger-caused by X).
 * @param nb          Neighborhood structure (adjacency list) of the lattice, defining the relationships between grid cells.
 * @param block       Block IDs for spatial block bootstrap (same length as x/y), indicating spatially contiguous subregions.
 * @param k           The number of neighbors to consider for symbolization (typically 3 or 5).
 * @param threads     Number of threads to use for parallel computation to speed up the bootstrap iterations.
 * @param boot        Number of bootstrap replicates to perform (default is 399).
 * @param base        Base of the logarithm for temporal lag scaling (default is 2, which corresponds to bits).
 * @param seed        Random seed for reproducibility of the results (default is 42).
 * @param symbolize   Flag indicating whether symbolization of the input grids is applied (default: true).
 * @param progressbar Whether to display a progress bar during the computation (default is true).
 *
 * @return A vector of four values:
 *         - scx: Observed statistic for causality from X → Y.
 *         - p_xy: Bootstrap-based p-value for X → Y.
 *         - scy: Observed statistic for causality from Y → X.
 *         - p_yx: Bootstrap-based p-value for Y → X.
 */
std::vector<double> SCT4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
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
    sc_bootstraps[n] = SCTSingle4Lattice(x_boot,y_boot,nb,static_cast<size_t>(std::abs(k)),base,symbolize);
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
  std::vector<double> sc = SCTSingle4Lattice(x,y,nb,static_cast<size_t>(std::abs(k)),base,symbolize);
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

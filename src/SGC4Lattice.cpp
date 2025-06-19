#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes directional spatial granger causality between two spatial variables
 * on a spatial lattice using spatial neighbor embeddings and quantized entropy measures.
 *
 * This function quantifies the asymmetric spatial granger causality between two
 * spatial variables `x` and `y`, both defined over a spatial lattice structure.
 * It adopts an information-theoretic framework based on symbolic (or continuous)
 * entropy estimation, incorporating spatial embedding through neighboring structures.
 *
 * Method Overview:
 * 1. Lattice-based Embedding:
 *    - For each spatial unit, embedding vectors `wx` and `wy` are generated using
 *      one-level spatial neighbors defined by the `nb` adjacency list.
 *
 * 2. Symbolization (Optional):
 *    - If `symbolize = true`, the inputs (`x`, `y`, `wx`, `wy`) are discretized into
 *      `k` symbolic categories prior to entropy computation.
 *      Otherwise, entropy is estimated directly from continuous values using kernel methods.
 *
 * 3. Entropy Computation:
 *    - The function calculates marginal and joint entropies including:
 *      H(x, wx), H(y, wy), H(wx), H(wy), H(wx, wy), H(wx, wy, x), H(wx, wy, y).
 *
 * 4. Directional Causality Strengths:
 *    - From x to y:
 *        sc_x_to_y = ((H(y, wy) - H(wy)) - (H(wx, wy, y) - H(wx, wy)))
 *    - From y to x:
 *        sc_y_to_x = ((H(x, wx) - H(wx)) - (H(wx, wy, x) - H(wx, wy)))
 *
 * 5. Normalization (Optional):
 *    - If `normalize = true`, the raw causality values are scaled by their respective
 *      baseline entropy gains to fall within the range [-1, 1]. This normalization
 *      enhances interpretability and comparability across different variable pairs
 *      or spatial datasets.
 *
 * Parameters:
 * - x: Input spatial variable `x` (vector of doubles).
 * - y: Input spatial variable `y` (same size as `x`).
 * - nb: Neighborhood list defining spatial adjacency (e.g., rook or queen contiguity).
 * - lib: A vector of indices representing valid neighbors to consider for each spatial unit.
 * - pred: A vector of indices specifying which elements to compute the spatial Granger causality.
 * - k: Number of discrete bins used for symbolization or KDE estimation.
 * - base: Logarithm base for entropy (default = 2, for bits).
 * - symbolize: Whether to apply symbolication for symbolic entropy (default = true).
 * - normalize: Whether to normalize causality values to the range [-1, 1] (default = false).
 *
 * Returns:
 *   A `std::vector<double>` of size 2:
 *     - [0] Estimated spatial granger causality from x to y
 *     - [1] Estimated spatial granger causality from y to x
 *   If `normalize = true`, both values are scaled to the range [-1, 1].
 */
std::vector<double> SGCSingle4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    size_t k,
    double base = 2,
    bool symbolize = true,
    bool normalize = false
){
  std::vector<double> wx;
  std::vector<std::vector<double>> Ex = GenLatticeEmbeddings(x,nb,1,1);
  for (const auto& row : Ex) {
    wx.insert(wx.end(), row.begin(), row.end());
  }

  std::vector<double> wy;
  std::vector<std::vector<double>> Ey = GenLatticeEmbeddings(y,nb,1,1);
  for (const auto& row : Ey) {
    wy.insert(wy.end(), row.begin(), row.end());
  }

  double Hwx, Hwy, Hxwx, Hywy, Hwxwy, Hwxwyx, Hwxwyy;

  if (symbolize){
    std::vector<double> sx = GenLatticeSymbolization(x,nb,lib,pred,k);
    std::vector<double> sy = GenLatticeSymbolization(y,nb,lib,pred,k);
    std::vector<double> swx = GenLatticeSymbolization(wx,nb,lib,pred,k);
    std::vector<double> swy = GenLatticeSymbolization(wy,nb,lib,pred,k);

    std::vector<std::vector<double>> sp_series(pred.size(),std::vector<double>(4));
    for (size_t i = 0; i < pred.size(); ++i){
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
    std::vector<std::vector<double>> sp_series(pred.size(),std::vector<double>(4));
    std::vector<double> wx_series(pred.size());
    std::vector<double> wy_series(pred.size());
    for (size_t i = 0; i < pred.size(); ++i){
      sp_series[i] = {x[pred[i]],y[pred[i]],wx[pred[i]],wy[pred[i]]}; // 0:x 1:y 2:wx 3:wy
      wx_series[i] = wx[pred[i]];
      wy_series[i] = wy[pred[i]];
    }

    Hxwx = CppJoinEntropy_Cont(sp_series, {0,2}, k, base, false); // H(x,wx)
    Hywy = CppJoinEntropy_Cont(sp_series, {1,3}, k, base, false); // H(y,wy)
    Hwx = CppEntropy_Cont(wx_series, k, base, false); // H(wx)
    Hwy = CppEntropy_Cont(wy_series, k, base, false); // H(wy)
    Hwxwy = CppJoinEntropy_Cont(sp_series,{2,3}, k, base, false); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Cont(sp_series,{0,2,3}, k, base, false); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Cont(sp_series,{1,2,3}, k, base, false); // H(wx,wy,y)
  }

  double sc_x_to_y,sc_y_to_x;
  if (normalize){
    // transformed to fall within [-1, 1]
    sc_x_to_y = ((Hywy - Hwy) - (Hwxwyy - Hwxwy)) / (Hywy - Hwy);
    sc_y_to_x = ((Hxwx - Hwx) - (Hwxwyx - Hwxwy)) / ((Hxwx - Hwx));
  } else {
    sc_x_to_y = (Hywy - Hwy) - (Hwxwyy - Hwxwy);
    sc_y_to_x = (Hxwx - Hwx) - (Hwxwyx - Hwxwy);
  }

  return {sc_x_to_y,sc_y_to_x};
}

/**
 * @brief Compute spatial granger causality for lattice data using spatial block bootstrap.
 *
 * This function estimates the directional spatial granger causality between two lattice variables `x` and `y`,
 * by applying a symbolic entropy-based method, and assesses the statistical significance of the causality using
 * spatial block bootstrap techniques. It calculates the causality in both directions: X → Y and Y → X.
 * Additionally, the function evaluates the significance of the estimated causality statistics by comparing them
 * to bootstrap realizations of the causality.
 *
 * The method involves the following steps:
 * - **Computation of true causality**: The function first calculates the spatial Granger causality statistic
 *   using the original lattice data `x` and `y`.
 * - **Spatial block bootstrap resampling**: The lattice values are resampled with spatial block bootstrapping.
 *   Each resample preserves local spatial structure and generates new bootstrap realizations of the causality statistic.
 * - **Estimation of causality for bootstrapped samples**: The causality statistic is estimated for each of the
 *   bootstrapped realizations, which involves calculating the symbolic entropy measures and their differences.
 * - **Empirical p-values**: The final p-values for both directional causality estimates (X → Y and Y → X) are
 *   derived by comparing the bootstrapped statistics with the true causality statistics.
 *
 * This approach accounts for spatial autocorrelation and allows the use of parallel processing for faster
 * bootstrap estimation. The spatial bootstrap method involves reshuffling lattice cells into spatial blocks,
 * preserving local dependencies, and calculating causality for each realization.
 *
 * @param x           Input vector for spatial variable x.
 * @param y           Input vector for spatial variable y (same length as x).
 * @param nb          Neighborhood list (e.g., queen or rook adjacency), used for embedding.
 * @param lib         A vector of indices representing valid neighbors to consider for each spatial unit.
 * @param pred        A vector of indices specifying which elements to compute the spatial Granger causality.
 * @param block       Vector indicating block assignments for spatial block bootstrapping.
 * @param k           Number of discrete bins used for symbolization or KDE estimation.
 * @param threads     Number of threads to use for parallel bootstrapping.
 * @param boot        Number of bootstrap iterations (default: 399).
 * @param base        Logarithmic base for entropy (default: 2, i.e., bits).
 * @param seed        Random seed for reproducibility (default: 42).
 * @param symbolize   Whether to apply symbolization before entropy computation (default: true).
 * @param normalize   Whether to normalize entropy values (optional, default: false).
 * @param progressbar Whether to display a progress bar during bootstrapping (default: true).
 *
 * @return A vector of four values:
 *         - sc_x_to_y: Estimated spatial granger causality from x to y.
 *         - p_x_to_y: Empirical p-value for x → y based on bootstrap distribution.
 *         - sc_y_to_x: Estimated spatial granger causality from y to x.
 *         - p_y_to_x: Empirical p-value for y → x based on bootstrap distribution.
 */
std::vector<double> SGC4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    const std::vector<int>& block,
    int k,
    int threads,
    int boot = 399,
    double base = 2,
    unsigned int seed = 42,
    bool symbolize = true,
    bool normalize = false,
    bool progressbar = true
){
  std::vector<std::vector<double>> sc_bootstraps(boot);

  // // Previous implementation may cause spurious correlations during randomization process
  // auto monte_boots = [&](int n){
  //   // Use different seed for each iteration to ensure different random samples
  //   unsigned int current_seed = seed + n;
  //   // Generate a spatial block bootstrap resample of indices
  //   std::vector<int> boot_indice = SpatialBlockBootstrap(block,current_seed);
  //   // Obtain the bootstrapped realization series
  //   std::vector<double> x_boot(x.size());
  //   std::vector<double> y_boot(y.size());
  //   for (size_t i = 0; i < boot_indice.size(); ++i){
  //     x_boot[i] = x[boot_indice[i]];
  //     y_boot[i] = y[boot_indice[i]];
  //   }
  //   // Estimate the bootstrapped realization of the spatial granger causality statistic
  //   sc_bootstraps[n] = SGCSingle4Lattice(x_boot,y_boot,nb,lib,pred,static_cast<size_t>(std::abs(k)),base,symbolize,normalize);
  // };

  // Prebuild RNG pool with seed sequence
  std::vector<std::mt19937> rng_pool(boot);
  for (int i = 0; i < boot; ++i) {
    std::seed_seq seq{static_cast<uint32_t>(seed), static_cast<uint32_t>(i)};
    rng_pool[i] = std::mt19937(seq);
  }

  auto monte_boots = [&](int n){
    // Use prebuilt rng instance
    std::vector<int> boot_indice = SpatialBlockBootstrapRNG(block, rng_pool[n]);

    std::vector<double> x_boot(x.size()), y_boot(y.size());
    for (size_t i = 0; i < boot_indice.size(); ++i){
      x_boot[i] = x[boot_indice[i]];
      y_boot[i] = y[boot_indice[i]];
    }

    sc_bootstraps[n] = SGCSingle4Lattice(
      x_boot, y_boot, nb, lib, pred,
      static_cast<size_t>(std::abs(k)), base,
      symbolize, normalize
    );
  };

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Parallel execution with progress bar if enabled
  if (progressbar) {
    RcppThread::ProgressBar bar(boot, 1);
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, boot, monte_boots, threads_sizet);
  }

  // Compute the original statistic (non-bootstrapped)
  std::vector<double> sc = SGCSingle4Lattice(
    x, y, nb, lib, pred,
    static_cast<size_t>(std::abs(k)), base,
    symbolize, normalize
  );

  double scx = sc[0], scy = sc[1];
  double b_xy = 0, b_yx = 0;
  for (size_t i = 0; i < sc_bootstraps.size(); ++i){
    if (sc_bootstraps[i][0] > scx) ++b_xy;
    if (sc_bootstraps[i][1] > scy) ++b_yx;
  }

  return {scx, b_xy / boot, scy, b_yx / boot};
}

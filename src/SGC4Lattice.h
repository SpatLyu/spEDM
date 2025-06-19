#ifndef SGC4Lattice_H
#define SGC4Lattice_H

#include <vector>
#include <cstdint>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

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
);

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
);

#endif // SGC4Lattice_H

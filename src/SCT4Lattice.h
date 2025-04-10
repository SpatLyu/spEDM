#ifndef SCT4Lattice_H
#define SCT4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

/*
 * @brief Computes directional spatial Granger causality between two spatial variables
 * on a spatial lattice data using spatial neighbor embeddings and quantized entropy measures.
 *
 * This function quantifies the asymmetric spatial Granger causality strength between two
 * spatial variables `x` and `y`, both defined over a spatial lattice data. It applies
 * a symbolic information-theoretic framework, incorporating spatial embedding via neighbor
 * structures and optional symbolization into discrete categories.
 *
 * Method Overview:
 * 1. Lattice-based Embedding:
 *    - For each spatial unit, embedded vectors `wx` and `wy` are constructed using 1-level
 *      neighbors based on the provided neighborhood list `nb`.
 *
 * 2. Symbolization (Optional):
 *    - If `symbolize = true`, the input vectors `x`, `y`, `wx`, and `wy` are discretized
 *      into `k` categories using lattice-based symbolization before entropy computation.
 *      Otherwise, continuous-valued entropy estimators are used.
 *
 * 3. Entropy Calculations:
 *    - The function computes marginal and joint entropies of various combinations of variables
 *      to evaluate symbolic causality. Specifically:
 *      - H(x, wx), H(y, wy), H(wx), H(wy), H(wx, wy), H(wx, wy, x), H(wx, wy, y)
 *
 * 4. Directional Causality Strengths:
 *    - From x to y:
 *      sc_x_to_y = (H(y, wy) - H(wy)) - (H(wx, wy, y) - H(wx, wy))
 *    - From y to x:
 *      sc_y_to_x = (H(x, wx) - H(wx)) - (H(wx, wy, x) - H(wx, wy))
 *    These values reflect the reduction in uncertainty in `y` (or `x`) when considering `x`
 *    (or `y`) and its spatial context.
 *
 * Parameters:
 * - x: Input spatial variable `x` (vector of doubles).
 * - y: Input spatial variable `y` (same size as `x`).
 * - nb: Neighborhood list defining spatial adjacency (e.g., rook or queen contiguity).
 * - k: Number of discrete bins used for symbolization or KDE estimation.
 * - base: Base of the logarithm for entropy (default = 2, for bits).
 * - symbolize: Whether to apply discretization for symbolic entropy (default = true).
 *
 * Returns:
 *    A `std::vector<double>` of size 2:
 *        - [0] Spatial Granger causality from `X` to `Y`.
 *        - [1] Spatial Granger causality from `Y` to `X`.
 */
std::vector<double> SCTSingle4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    double base = 2,
    bool symbolize = true
);

/**
 * @brief Compute spatial Granger causality for lattice data using spatial block bootstrap.
 *
 * This function estimates the directional spatial Granger causality between two lattice variables `x` and `y`,
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
 * @param block       Vector indicating block assignments for spatial block bootstrapping.
 * @param k           Number of discrete bins used for symbolization or KDE estimation.
 * @param threads     Number of threads to use for parallel bootstrapping.
 * @param boot        Number of bootstrap iterations (default: 399).
 * @param base        Logarithmic base for entropy (default: 2, i.e., bits).
 * @param seed        Random seed for reproducibility (default: 42).
 * @param symbolize   Whether to apply symbolization before entropy computation (default: true).
 * @param progressbar Whether to display a progress bar during bootstrapping (default: true).
 *
 * ## Returns:
 * A std::vector<double> of length 4:
 * - [0] Observed symbolic causality from x to y (sc_x_to_y).
 * - [1] Empirical p-value for x → y based on bootstrap distribution.
 * - [2] Observed symbolic causality from y to x (sc_y_to_x).
 * - [3] Empirical p-value for y → x based on bootstrap distribution.
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
);

#endif // SCT4Lattice_H

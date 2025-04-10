#ifndef SCT4Grid_H
#define SCT4Grid_H

#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

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
);

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
);

#endif // SCT4Grid_H

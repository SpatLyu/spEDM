#ifndef SCT4Grid_H
#define SCT4Grid_H

#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

/**
 * @brief Compute the spatial Granger causality statistic between two grid-structured variables.
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
 * @param x     2D grid of variable X (as a vector of rows).
 * @param y     2D grid of variable Y (same size as x).
 * @param k     Neighborhood window size for local symbolization.
 * @param base  Logarithm base used in entropy computation (default is 2 for bits).
 *
 * @return A vector with two values:
 *         - sc_x_to_y: Spatial Granger causality from X to Y.
 *         - sc_y_to_x: Spatial Granger causality from Y to X.
 */
std::vector<double> SCTSingle4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    size_t k,
    double base = 2
);

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
);

#endif // SCT4Grid_H

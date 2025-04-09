#ifndef SCT4Grid_H
#define SCT4Grid_H

#include <vector>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

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
 *
 * @return A std::vector<double> of two values:
 *         - sc_x_to_y: Directional spatial Granger causality from X to Y.
 *         - sc_y_to_x: Directional spatial Granger causality from Y to X.
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
);

#endif // SCT4Grid_H

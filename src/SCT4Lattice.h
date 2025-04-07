#ifndef SCT4Lattice_H
#define SCT4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

/**
 * @brief Computes the directional symbolic causality strength between two spatial variables
 *        based on lattice neighborhood structure and optional symbolization.
 *
 * This function implements a symbolic causality test for spatial data over lattice structures.
 * It evaluates whether variable x symbolically Granger-causes variable y using entropy-based
 * measures on either the raw or symbolized values of x and y, depending on the `symbolize` flag.
 *
 * Specifically, the function calculates:
 *   - H(y): the entropy of y (or symbolized y if `symbolize` is true).
 *   - H(y | x): the conditional entropy of y given x (or symbolized x if `symbolize` is true).
 *   - H(x): the entropy of x (or symbolized x).
 *   - H(x | y): the conditional entropy of x given y (or symbolized y).
 * Then it computes:
 *   - sc_x_to_y = H(y) - H(y | x): symbolic causality strength from x to y.
 *   - sc_y_to_x = H(x) - H(x | y): symbolic causality strength from y to x.
 *
 * @param x         Input vector of variable x, aligned with a lattice/grid.
 * @param y         Input vector of variable y, aligned with the same lattice/grid.
 * @param nb        Neighborhood list for each lattice location (e.g., rook or queen neighbors).
 * @param k         Number of neighbors used in lattice symbolization.
 * @param base      Logarithm base for entropy calculation (default is 2).
 * @param symbolize Whether to perform lattice-based symbolization on x and y before entropy computation.
 *
 * @return A std::vector<double> of length 2:
 *         [0] sc_x_to_y — symbolic causality strength from x to y;
 *         [1] sc_y_to_x — symbolic causality strength from y to x.
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
    bool progressbar = true
);

#endif // SCT4Lattice_H

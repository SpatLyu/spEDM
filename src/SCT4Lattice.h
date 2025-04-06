#ifndef SCT4Lattice_H
#define SCT4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

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
);

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
);

#endif // SCT4Lattice_H

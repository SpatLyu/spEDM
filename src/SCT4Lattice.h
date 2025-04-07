#ifndef SCT4Lattice_H
#define SCT4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

/**
 * @brief Computes directional symbolic causality strength between two spatial variables
 *        over a lattice using neighbor-based embedding and mandatory symbolization.
 *
 * This function evaluates the symbolic directional influence (causality) between two spatial variables
 * `x` and `y` that are aligned on a regular lattice or grid. It uses a symbolic information-theoretic
 * approach, incorporating both spatial embedding and discrete symbolization based on neighborhood
 * structure.
 *
 * ## Method Overview:
 * 1. **Lattice Embedding**:
 *    - For each location, generate embedded (lagged) spatial vectors `wx` and `wy` using the neighborhood list `nb`.
 *    - The embeddings are based on a window of size 1 in space (i.e., 1-level neighbors).
 *
 * 2. **Lattice Symbolization**:
 *    - Symbolize the original variables `x`, `y` and their embeddings `wx`, `wy` using lattice-based quantization
 *      into `k` discrete categories.
 *
 * 3. **Entropy Computation**:
 *    - Calculate marginal and joint entropies using the symbolized variables:
 *      - \( H(x, wx), H(y, wy), H(wx), H(wy), H(wx, wy), H(wx, wy, x), H(wx, wy, y) \)
 *
 * 4. **Symbolic Causality Strength**:
 *    - From x to y:
 *      \[
 *      sc_{x \rightarrow y} = [H(y, wy) - H(wy)] - [H(wx, wy, y) - H(wx, wy)]
 *      \]
 *    - From y to x:
 *      \[
 *      sc_{y \rightarrow x} = [H(x, wx) - H(wx)] - [H(wx, wy, x) - H(wx, wy)]
 *      \]
 *    These quantities represent the reduction in uncertainty of `y` (or `x`) due to `x` (or `y`)
 *    in a symbolic and spatial context.
 *
 * @param x     Input vector of spatial variable x (must match the lattice structure).
 * @param y     Input vector of spatial variable y (same size as x).
 * @param nb    Neighborhood structure for each location (e.g., rook or queen adjacency).
 * @param k     Number of discrete bins used in symbolization of each variable.
 * @param base  Base of the logarithm used for entropy computation (default is 2, i.e., bits).
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
    double base = 2
);

/**
 * @brief Computes symbolic causality strength on a spatial lattice using bootstrap for significance testing.
 *
 * This function evaluates the directional symbolic causality strength between two spatial variables `x` and `y`,
 * organized on a regular lattice structure with known neighborhood relations (`nb`). It uses symbolic entropy-based
 * causality measures, extended to the spatial domain using lattice embedding and neighborhood-driven symbolization.
 *
 * To assess the statistical significance of the observed causality, a spatial block bootstrap approach is employed.
 * This method resamples spatial blocks (as specified in `block`) to preserve spatial dependence structures in the
 * data, and computes bootstrap replicates of the causality statistic. Empirical p-values are derived by comparing
 * observed statistics against this bootstrap distribution.
 *
 * Core steps:
 * 1. **Causality Computation**: The directional symbolic causality strength from `x` to `y` and from `y` to `x` is
 *    computed using the `SCTSingle4Lattice()` function, which applies lattice embedding, (optional) symbolization,
 *    and entropy-based measures to evaluate uncertainty reduction.
 * 2. **Spatial Block Bootstrap**: For each of `boot` replicates, spatially contiguous block indices (based on `block`)
 *    are resampled to create bootstrap versions of `x` and `y`. Causality statistics are recomputed for each.
 * 3. **Significance Evaluation**: p-values are computed by counting how often the bootstrapped statistics exceed
 *    the observed statistic in each direction.
 *
 * @param x           Input spatial variable `x`, aligned with lattice cells.
 * @param y           Input spatial variable `y`, aligned with lattice cells.
 * @param nb          Neighborhood list, specifying for each cell its adjacent cells.
 * @param block       Vector of block IDs for spatial bootstrap (one per cell), indicating spatial clusters.
 * @param k           Number of bins used for symbolization (typically 3 to 5).
 * @param threads     Number of threads to use for parallel bootstrap computation.
 * @param boot        Number of bootstrap replicates (default: 399).
 * @param base        Logarithmic base used in entropy calculation (default: 2).
 * @param seed        Random seed to ensure reproducibility (default: 42).
 * @param progressbar Whether to display a progress bar for bootstrap iterations (default: true).
 *
 * @return A std::vector<double> of size 4:
 *         - [0] Observed causality strength from `x` to `y`.
 *         - [1] Bootstrap-based p-value for `x` → `y`.
 *         - [2] Observed causality strength from `y` to `x`.
 *         - [3] Bootstrap-based p-value for `y` → `x`.
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

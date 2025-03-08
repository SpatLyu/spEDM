#ifndef GCCM4Lattice_H
#define GCCM4Lattice_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

/*
 * Perform GCCM on a single lib and pred for lattice data.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: Spatial cross-section series used as the target (should align with x_vectors).
 *   - lib_size: Size of the library used for cross mapping.
 *   - max_lib_size: Maximum size of the library.
 *   - possible_lib_indices: Indices of possible library states.
 *   - pred_indices: A boolean vector indicating which states to use for prediction.
 *   - b: Number of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *   - threads: The number of threads to use for parallel processing.
 *
 * Returns:
 *   A vector of pairs, where each pair consists of:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho) between predicted and actual values.
 */
std::vector<std::pair<int, double>> GCCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,
    const std::vector<double>& y,
    int lib_size,
    int max_lib_size,
    const std::vector<int>& possible_lib_indices,
    const std::vector<bool>& pred_indices,
    int b,
    bool simplex,
    double theta,
    size_t threads
);

/**
 * Performs GCCM on a spatial lattice data.
 *
 * Parameters:
 * - x: Spatial cross-section series used as the predict variable (**cross mapping from**).
 * - y: Spatial cross-section series used as the target variable (**cross mapping to**).
 * - nb_vec: A nested vector containing neighborhood information for lattice data.
 * - lib_sizes: A vector specifying different library sizes for GCCM analysis.
 * - lib: A vector specifying the library indices (1-based in R, converted to 0-based in C++).
 * - pred: A vector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 * - E: Embedding dimension for attractor reconstruction.
 * - tau: the step of spatial lags for prediction.
 * - b: Number of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - progressbar: Boolean flag to indicate whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean cross-mapping correlation.
 *      - The statistical significance of the correlation.
 *      - The lower bound of the confidence interval.
 *      - The upper bound of the confidence interval.
 */
std::vector<std::vector<double>> GCCM4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb_vec,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    bool progressbar
);

#endif // GCCM4Lattice_H

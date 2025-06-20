#ifndef IntersectionCardinality_H
#define IntersectionCardinality_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include "CppStats.h"
#include "spEDMDataStruct.h"
#include <RcppThread.h>

/**
 * Computes intersection-based mapping ratio sequences between two neighbor graphs
 * for use in Cross Mapping Cardinality (CMC) or similar causal inference frameworks.
 *
 * Parameters:
 *   neighborsX     - Precomputed sorted neighbor indices for embedding X
 *   neighborsY     - Precomputed sorted neighbor indices for embedding Y
 *   lib_size       - Size of the moving library used in mapping
 *   lib_indices    - Global indices from which to draw the sliding libraries
 *   pred_indices   - Indices at which to perform prediction (evaluation points)
 *   num_neighbors  - Number of neighbors used for mapping (after exclusion)
 *   n_excluded     - Number of nearest neighbors to exclude from the front
 *   threads        - Number of parallel threads for computation
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of IntersectionRes structures, each containing the average intersection
 *   ratio sequence (IC curve) for a different starting point of the moving library.
 *   If lib_size == lib_indices.size(), returns a single result using full library.
 *
 * Notes:
 *   - Neighbor lists must use std::numeric_limits<size_t>::max() to indicate invalid entries.
 *   - This function assumes that the neighbor vectors are sorted by ascending distance.
 *   - Use in combination with AUC computation to assess causal strength.
 */
std::vector<IntersectionRes> IntersectionCardinalitySingle(
    const std::vector<std::vector<size_t>>& neighborsX,
    const std::vector<std::vector<size_t>>& neighborsY,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    size_t num_neighbors,
    size_t n_excluded,
    size_t threads,
    int parallel_level
);

/*
 * Computes the Intersection Cardinality (IC) scores
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   lib: Library index vector (1-based in R, converted to 0-based).
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Number of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Number of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   - A vector representing the intersection cardinality (IC) scores, normalized between [0, 1].
 */
std::vector<double> IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int num_neighbors,
    int n_excluded,
    int threads,
    bool progressbar);

#endif // IntersectionCardinality_H

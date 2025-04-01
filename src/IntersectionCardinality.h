#ifndef IntersectionCardinality_H
#define IntersectionCardinality_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <unordered_set>
#include "CppStats.h"
#include <RcppThread.h>

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

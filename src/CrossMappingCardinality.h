#ifndef CrossMappingCardinality_H
#define CrossMappingCardinality_H

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
 * Computes the Intersection Cardinality (IC) causal strength score.
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Number of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Number of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   - A double representing the IC causal strength score, normalized between [0,1].
 */
double IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    int num_neighbors,
    int n_excluded,
    int threads,
    bool progressbar);

/**
 * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   A vector of normalized CMC causal strength scores in the range [0,1], corresponding to each value in num_neighbors.
 */
std::vector<double> CrossMappingCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    const std::vector<int>& num_neighbors,
    const std::vector<int>& n_excluded,
    int threads,
    bool progressbar);

/**
 * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   A vector of normalized CMC causal strength scores in the range [0,1], corresponding to each value in num_neighbors.
 */
std::vector<double> CrossMappingCardinality2(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    const std::vector<int>& num_neighbors,
    const std::vector<int>& n_excluded,
    int threads,
    bool progressbar);

#endif // CrossMappingCardinality_H

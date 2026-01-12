#ifndef MultiViewEmbedding_H
#define MultiViewEmbedding_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "CppStats.h"
#include "SimplexProjection.h"
// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.

/**
 * Computes the multi-view embedding by evaluating multiple feature embeddings using simplex projection,
 * selecting top-performing embeddings, and aggregating their contributions.
 *
 * Parameters:
 *   - vectors: 2D vector where each row represents a sample and each column a feature.
 *   - target: Target spatial cross sectional series aligned with the samples in vectors.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - num_neighbors: Number of neighbors used for simplex projection.
 *   - top_num: Number of top-performing reconstructions to select.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A vector<double> where each element is the predict value from selected embeddings columns.
 */
std::vector<double> MultiViewEmbedding(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    int top_num = 3,
    int dist_metric = 2,
    int dist_average = true,
    int threads = 8
);

/**
 * @brief Computes the multi-view embedding from stacked (3D) feature embeddings using simplex projection.
 *
 * This overload supports a 3D nested vector structure where each outer element
 * represents a feature stack or multivariate embedding set. Each stack contains
 * a 2D matrix (samples × embedded dimensions). The function evaluates each stack
 * independently using simplex projection (via SimplexBehavior), ranks the stacks
 * by their predictive performance, and constructs a new 3D subset of the top
 * performing stacks for final prediction.
 *
 * Unlike the 2D overload, this version preserves the structural independence of
 * each selected embedding stack, forming a 3D tensor of selected embeddings that
 * is passed to SimplexProjectionPrediction for final multi-view prediction.
 *
 * @param vectors       3D vector: [stack][sample][dimension], representing multiple feature stacks.
 * @param target        Target spatial cross-sectional series aligned with sample indices.
 * @param lib_indices   Indices for the library (training) set.
 * @param pred_indices  Indices for the prediction set.
 * @param num_neighbors Number of neighbors used for simplex projection.
 * @param top_num       Number of top-performing stacks to select.
 * @param dist_metric   Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average  Whether to average distances by the number of valid vector components.
 * @param threads       Number of threads used for parallel computation.
 *
 * @return A vector<double> containing predicted values derived from the selected embedding stacks.
 */
std::vector<double> MultiViewEmbedding(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    int top_num = 3,
    int dist_metric = 2,
    int dist_average = true,
    int threads = 8
);

#endif // MultiViewEmbedding_H

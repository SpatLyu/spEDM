#ifndef MultiViewEmbedding_H
#define MultiViewEmbedding_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "CppStats.h"
#include "SimplexProjection.h"
#include <RcppThread.h>

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
    int num_neighbors,
    int top_num,
    int threads
);

#endif // MultiViewEmbedding_H

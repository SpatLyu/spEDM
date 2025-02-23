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
 * Computes the multi-view embedding by evaluating multiple feature subsets using simplex projection,
 * selecting top-performing subsets, and aggregating their contributions.
 *
 * Parameters:
 *   - vectors: 2D vector where each row represents a sample and each column a feature.
 *   - target: Target spatial cross sectional series aligned with the samples in vectors.
 *   - lib_indices: Boolean flag indicating which samples to use for neighbor search.
 *   - pred_indices: Boolean flag indicating which samples to predict.
 *   - num_neighbors: Number of neighbors used for simplex projection.
 *   - top_num: Number of top-performing reconstructions to select.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A vector<double> where each element is the aggregated average of predict value from selected feature columns.
 */
std::vector<double> MultiViewEmbedding(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors,
    int top_num,
    int threads
);

#endif // MultiViewEmbedding_H

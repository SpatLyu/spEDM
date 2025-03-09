#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/*
 * Computes predictions using the simplex projection method based on state-space reconstruction.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Spatial cross sectional series used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of T/F values indicating which states to predict from.
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A vector<double> containing predicted target values (target_pred).
 */
std::vector<double> SimplexProjectionPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
) {
  // Convert num_neighbors to size_t
  size_t num_neighbors_sizet = static_cast<size_t>(num_neighbors);

  // Setup output
  std::vector<double> pred(target.size(), std::numeric_limits<double>::quiet_NaN());

  // no neighbor to use, return all nan
  if (num_neighbors <= 0){
    return pred;
  }

  // Count the number of true values in lib_indices
  size_t lib_count = std::count(lib_indices.begin(), lib_indices.end(), true);

  // no library to use, return all nan
  if (lib_count == 0){
    return pred;
  }

  // // If the number of true values is less than num_neighbors, return NaN vector
  // if (lib_count < num_neighbors_sizet) {
  //   return pred;
  // }

  // Make predictions
  for (size_t p = 0; p < pred_indices.size(); ++p) {
    if (!pred_indices[p]) continue;

    // Create a local copy of lib_indices to modify
    std::vector<bool> local_lib_indices = lib_indices;
    bool temp_lib = local_lib_indices[p];
    local_lib_indices[p] = false;
    std::vector<size_t> libs;
    for (size_t i = 0; i < local_lib_indices.size(); ++i) {
      if (local_lib_indices[i]) libs.push_back(i);
    }

    // Handle the case where libs is empty
    if (libs.empty()) {
      pred[p] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    // Adjust num_neighbors_sizet if it exceeds libs.size()
    if (num_neighbors_sizet > libs.size()) {
      num_neighbors_sizet = libs.size();
    }

    // Compute distances
    std::vector<double> distances;
    for (size_t i : libs) {
      double sum_sq = 0.0;
      double sum_na = 0.0;
      for (size_t j = 0; j < vectors[p].size(); ++j) {
        if (!std::isnan(vectors[i][j]) && !std::isnan(vectors[p][j])) {
          sum_sq += std::pow(vectors[i][j] - vectors[p][j], 2);
          sum_na += 1.0;
        }
      }

      if (sum_na > 0) {
        distances.push_back(std::sqrt(sum_sq / sum_na));
      } else {
        distances.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }

    // Find nearest neighbors
    std::vector<size_t> neighbors(distances.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);
    std::partial_sort(neighbors.begin(), neighbors.begin() + num_neighbors_sizet, neighbors.end(),
                      [&](size_t a, size_t b) {
                        return (distances[a] < distances[b]) ||
                               (distances[a] == distances[b] && a < b);
                        });

    double min_distance = distances[neighbors[0]];

    // Compute weights
    std::vector<double> weights(num_neighbors_sizet);
    if (min_distance == 0) { // Perfect match
      std::fill(weights.begin(), weights.end(), 0.000001);
      for (size_t i = 0; i < num_neighbors_sizet; ++i) {
        if (distances[neighbors[i]] == 0) weights[i] = 1.0;
      }
    } else {
      for (size_t i = 0; i < num_neighbors_sizet; ++i) {
        weights[i] = std::exp(-distances[neighbors[i]] / min_distance);
        if (weights[i] < 0.000001) weights[i] = 0.000001;
      }
    }
    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    // Make prediction(use inner product or iterate element-wise for computation)

    // std::vector<double> target_neighbors;
    // for (size_t i = 0; i < num_neighbors_sizet; ++i) {
    //   target_neighbors.push_back(target[libs[neighbors[i]]]);
    // }
    // double prediction = std::inner_product(weights.begin(), weights.end(), target_neighbors.begin(), 0.0);

    double prediction = 0.0;
    for (size_t i = 0; i < num_neighbors_sizet; ++i) {
      prediction += weights[i] * target[libs[neighbors[i]]];
    }

    pred[p] = prediction / total_weight;

    // Restore the original lib_indices state
    local_lib_indices[p] = temp_lib;
  }

  return pred;
}

/*
 * Computes the Pearson correlation coefficient (rho) using the simplex projection prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Spatial cross sectional series used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of T/F values indicating which states to use for prediction.
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A double representing the Pearson correlation coefficient (rho) between the predicted and actual target values.
 */
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
) {
  double rho = std::numeric_limits<double>::quiet_NaN();

  // Call SimplexProjectionPrediction to get the prediction results
  std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    rho = PearsonCor(target_pred, target, true);
  }
  return rho;
}

/*
 * Computes the simplex projection and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Spatial cross sectional series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
 *   - pred_indices: Vector of T/F values (which states to predict from).
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A vector<double> containing:
 *     - Pearson correlation coefficient (PearsonCor)
 *     - Mean absolute error (MAE)
 *     - Root mean squared error (RMSE)
 */
std::vector<double> SimplexBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
) {
  // Initialize PearsonCor, MAE, and RMSE
  double pearson = std::numeric_limits<double>::quiet_NaN();
  double mae = std::numeric_limits<double>::quiet_NaN();
  double rmse = std::numeric_limits<double>::quiet_NaN();

  // Call SimplexProjectionPrediction to get the prediction results
  std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    // Compute PearsonCor, MAE, and RMSE
    pearson = PearsonCor(target_pred, target, true);
    mae = CppMAE(target_pred, target, true);
    rmse = CppRMSE(target_pred, target, true);
  }

  // Return the three metrics as a vector
  return {pearson, mae, rmse};
}

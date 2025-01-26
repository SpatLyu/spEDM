#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

// Function to compute the 'simplex projection' prediction
// Returns a pair of vectors: (target_pred, pred_pred)
std::pair<std::vector<double>, std::vector<double>> SimplexProjectionPrediction(
    const std::vector<std::vector<double>>& vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<bool>& lib_indices,             // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,            // Vector of T/F values (which states to predict from)
    int num_neighbors                                 // Number of neighbors to use for simplex projection
) {
  // Convert num_neighbors to size_t
  size_t num_neighbors_sizet = static_cast<size_t>(num_neighbors);

  // Setup output
  std::vector<double> pred(target.size(), std::numeric_limits<double>::quiet_NaN());

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
      distances.push_back(std::sqrt(sum_sq / sum_na));
    }

    // Find nearest neighbors
    std::vector<size_t> neighbors(distances.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);
    // Check if num_neighbors_sizet exceeds the size of neighbors
    // If so, adjust num_neighbors_sizet to the maximum allowed value (neighbors.size())
    if (num_neighbors_sizet > neighbors.size()) {
      num_neighbors_sizet = neighbors.size();
    }
    std::partial_sort(neighbors.begin(), neighbors.begin() + num_neighbors_sizet, neighbors.end(),
                      [&](size_t a, size_t b) { return distances[a] < distances[b]; });

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

    // Make prediction
    double prediction = 0.0;
    for (size_t i = 0; i < num_neighbors_sizet; ++i) {
      prediction += weights[i] * target[libs[neighbors[i]]];
    }
    pred[p] = prediction / total_weight;

    // Restore the original lib_indices state
    local_lib_indices[p] = temp_lib;
  }

  // Extract the target and prediction values for the prediction indices
  std::vector<double> target_pred;
  std::vector<double> pred_pred;
  for (size_t i = 0; i < pred_indices.size(); ++i) {
    if (pred_indices[i]) {
      target_pred.push_back(target[i]);
      pred_pred.push_back(pred[i]);
    }
  }

  // Return the pair of vectors (target_pred, pred_pred)
  return std::make_pair(target_pred, pred_pred);
}

// Rho value by the 'simplex projection' prediction
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<bool>& lib_indices,             // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,            // Vector of T/F values (which states to predict from)
    int num_neighbors                                 // Number of neighbors to use for simplex projection
) {
  std::pair<std::vector<double>, std::vector<double>> result = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);
  std::vector<double> target_pred = result.first;
  std::vector<double> pred_pred = result.second;

  // Return the Pearson correlation coefficient
  return PearsonCor(target_pred, pred_pred, true);
}

// Description: Computes the simplex projection and returns a vector containing
//              Pearson correlation coefficient (PearsonCor), mean absolute error (MAE),
//              and root mean squared error (RMSE).
// Parameters:
//   - vectors: Reconstructed state-space (each row is a separate vector/state).
//   - target: Spatial cross-section series to be used as the target (should align with vectors).
//   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
//   - pred_indices: Vector of T/F values (which states to predict from).
//   - num_neighbors: Number of neighbors to use for simplex projection.
// Returns: A vector<double> containing {PearsonCor, MAE, RMSE}.
std::vector<double> SimplexBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
) {
  // Call SimplexProjectionPrediction to get the prediction results
  std::pair<std::vector<double>, std::vector<double>> result = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);
  std::vector<double> target_pred = result.first;
  std::vector<double> pred_pred = result.second;

  // Compute PearsonCor, CppMAE, and CppRMSE
  double pearson = PearsonCor(target_pred, pred_pred, true);
  double mae = CppMAE(target_pred, pred_pred);
  double rmse = CppRMSE(target_pred, pred_pred);

  // Return the three metrics as a vector
  return {pearson, mae, rmse};
}

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"
#include "CppUtils.h"

// Function to compute the simplex projection
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                // Time series to be used as the target (should line up with vectors)
    const std::vector<bool>& lib_indices,             // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,            // Vector of T/F values (which states to predict from)
    int num_neighbors                                // Number of neighbors to use for simplex projection
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

  // Return the Pearson correlation coefficient
  return PearsonCor(target_pred, pred_pred, true);
}

// Function to compute GCCMSingle4Lattice
std::vector<std::pair<int, double>> GCCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                      // Time series to be used as the target (should line up with vectors)
    const std::vector<bool>& lib_indices,              // Vector of T/F values (which states to include when searching for neighbors)
    int lib_size,                                      // Size of the library
    int max_lib_size,                                  // Maximum size of the library
    const std::vector<int>& possible_lib_indices,      // Indices of possible library states
    const std::vector<bool>& pred_indices,             // Vector of T/F values (which states to predict from)
    int b                                              // Number of neighbors to use for simplex projection
) {
  int n = x_vectors.size();
  std::vector<std::pair<int, double>> x_xmap_y;

  if (lib_size == max_lib_size) { // No possible library variation if using all vectors
    std::vector<bool> lib_indices(n, false);
    for (int idx : possible_lib_indices) {
      lib_indices[idx] = true;
    }

    // Run cross map and store results
    double rho = SimplexProjection(x_vectors, y, lib_indices, pred_indices, b);
    x_xmap_y.emplace_back(lib_size, rho);
  } else {
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<bool> lib_indices(n, false);
      // Setup changing library
      if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
        for (int i = start_lib; i < max_lib_size; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
      }

      // Run cross map and store results
      double rho = SimplexProjection(x_vectors, y, lib_indices, pred_indices, b);
      x_xmap_y.emplace_back(lib_size, rho);
    }
  }

  return x_xmap_y;
}

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "HelperFuns.h"
#include "CppStats.h"
#include "CppUtils.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

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

// Function to compute GCCMLattice
std::vector<std::vector<double>> GCCMLattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                      // Time series to cross map to
    const std::vector<int>& lib_sizes,                 // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& lib,       // Matrix (n x 2) using n sequences of data to construct libraries
    const std::vector<std::pair<int, int>>& pred,      // Matrix (n x 2) using n sequences of data to predict from
    int E,                                             // Number of dimensions for the attractor reconstruction
    int tau = 1,                                       // Time lag for the lagged-vector construction
    int b = 0                                          // Number of nearest neighbors to use for prediction
) {
  int n = x_vectors.size();
  b = E + 1; // Set b to E + 1 if not provided

  // Setup pred_indices
  std::vector<bool> pred_indices(n, false);
  for (const auto& p : pred) {
    int row_start = p.first + (E - 1) * tau;
    int row_end = p.second;
    if (row_end > row_start) {
      std::fill(pred_indices.begin() + row_start, pred_indices.begin() + row_end + 1, true);
    }
  }

  // Setup lib_indices
  std::vector<bool> lib_indices(n, false);
  for (const auto& l : lib) {
    int row_start = l.first + (E - 1) * tau;
    int row_end = l.second;
    if (row_end > row_start) {
      std::fill(lib_indices.begin() + row_start, lib_indices.begin() + row_end + 1, true);
    }
  }

  int max_lib_size = std::accumulate(lib_indices.begin(), lib_indices.end(), 0); // Maximum lib size
  std::vector<int> possible_lib_indices;
  for (int i = 0; i < n; ++i) {
    if (lib_indices[i]) {
      possible_lib_indices.push_back(i);
    }
  }

  // Make sure max lib size not exceeded and remove duplicates;
  // Ensure the minimum value in unique_lib_sizes is E + 2
  std::vector<int> unique_lib_sizes;
  std::unique_copy(lib_sizes.begin(), lib_sizes.end(), std::back_inserter(unique_lib_sizes),
                   [&](int a, int b) { return a == b; });
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::max(size, E + 2); });

  std::vector<std::pair<int, double>> x_xmap_y;

  // // Sequential version of the for loop
  // for (int lib_size : unique_lib_sizes) {
  //   auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b);
  //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
  // }

  // Perform the operations using RcppThread
  RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
  RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
    int lib_size = unique_lib_sizes[i];
    auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b);
    x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    bar++;
  });

  // Group by the first int and compute the mean
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  std::vector<std::vector<double>> final_results;
  for (const auto& group : grouped_results) {
    double mean_value = CppMean(group.second, true);
    final_results.push_back({static_cast<double>(group.first), mean_value});
  }

  // Calculate significance and confidence interval for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho = final_results[i][1];
    double significance = CppSignificance(rho, n);
    std::vector<double> confidence_interval = CppConfidence(rho, n);

    final_results[i].push_back(significance);
    final_results[i].push_back(confidence_interval[0]);
    final_results[i].push_back(confidence_interval[1]);
  }

  return final_results;
}

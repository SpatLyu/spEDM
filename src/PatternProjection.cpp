#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <utility>
#include <limits>
#include "DataStruct.h"

PatternProjectionRes PatternProjectionSingle(
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& Dx,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = std::numeric_limits<int>::min(),
    int zero_tolerance = std::numeric_limits<int>::min()
) {
  const size_t& n_row = SMy.size();
  const size_t& n_col = SMy[0].size();

  if (num_neighbors == std::numeric_limits<int>::min()) {
    num_neighbors = static_cast<int>(n_col + 2);
  }

  if (zero_tolerance == std::numeric_limits<int>::min()) {
    zero_tolerance = static_cast<int>(n_col);
  }

  std::vector<std::vector<double>> predSMy(
      n_row,std::vector<double>(n_col,std::numeric_limits<double>::quiet_NaN()));

  std::vector<std::vector<std::uint8_t>> predPMy(
      n_row,std::vector<std::uint8_t>(n_col, 0));

  for (size_t pi = 0; pi < pred_indices.size(); ++pi) {
    int p = pred_indices[pi];

    // Access distances only for valid vector pairs (exclude NaNs)
    std::vector<double> distances;
    distances.reserve(lib_indices.size());
    // keep track of libs corresponding to valid distances
    std::vector<int> valid_libs;
    valid_libs.reserve(lib_indices.size());

    for (size_t i : lib_indices) {
      if (!std::isnan(Dx[i][p])) {
        distances.push_back(Dx[i][p]);
        valid_libs.push_back(i);
      }
    }

    // If no valid distances found, prediction is NaN
    if (distances.empty()) {
      continue;
    }

    // Adjust number of neighbors to available valid libs
    size_t k = std::min(static_cast<size_t>(num_neighbors), distances.size());

    // Prepare indices for sorting
    std::vector<size_t> neighbors(distances.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);

    // Partial sort to find k nearest neighbors by distance
    std::partial_sort(
      neighbors.begin(), neighbors.begin() + k, neighbors.end(),
      [&](size_t a, size_t b) {
        return (distances[a] < distances[b]) ||
          (distances[a] == distances[b] && a < b);
      });

    double min_distance = distances[neighbors[0]];

    // Compute weights for neighbors
    std::vector<double> weights(k);
    if (min_distance == 0.0) {  // Perfect match found
      std::fill(weights.begin(), weights.end(), 0.000001);
      for (size_t i = 0; i < k; ++i) {
        if (distances[neighbors[i]] == 0.0) {
          weights[i] = 1.0;
        }
      }
    } else {
      for (size_t i = 0; i < k; ++i) {
        weights[i] = std::exp(-distances[neighbors[i]] / min_distance);
        if (weights[i] < 0.000001) {
          weights[i] = 0.000001;
        }
      }
    }

    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    // Calculate weighted prediction, ignoring any NaN targets
    // (No NaNs here, as NaN values in the corresponding components of lib and pred are excluded in advance.)
    double prediction = 0.0;
    for (size_t i = 0; i < k; ++i) {
      prediction += weights[i] * target[valid_libs[neighbors[i]]];
    }

    pred[p] = prediction / total_weight;
  }

}

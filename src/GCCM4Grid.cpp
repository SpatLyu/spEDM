#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <limits>
#include "CppStats.h"

// Function to locate the index in a 2D grid
int locate(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
}

// Function to calculate the distance between embeddings, handling NaN values
std::vector<double> DistCom(const std::vector<std::vector<std::vector<double>>>& embeddings,
                            const std::vector<int>& libs, int p) {
  std::vector<double> distances;
  distances.reserve(libs.size());

  // Iterate over each embedding
  for (const auto& embedding : embeddings) {
    // Get the p-th element of the current embedding
    const std::vector<double>& emd_p = embedding[p];

    // Calculate the absolute difference between the p-th element and each library element
    for (int lib : libs) {
      const std::vector<double>& emd_lib = embedding[lib];
      double distance = 0.0;
      bool has_nan = false; // Flag to check if any NaN is encountered

      // Compute the distance while checking for NaN
      for (size_t i = 0; i < emd_p.size(); ++i) {
        if (std::isnan(emd_p[i]) || std::isnan(emd_lib[i])) {
          has_nan = true;
          break;
        }
        distance += std::abs(emd_p[i] - emd_lib[i]);
      }

      // If NaN is encountered, skip this distance
      if (has_nan) {
        distances.push_back(std::numeric_limits<double>::quiet_NaN());
      } else {
        distances.push_back(distance);
      }
    }
  }

  // Calculate row-wise mean, excluding NaN values
  std::vector<double> row_means(libs.size(), 0.0);
  for (size_t i = 0; i < libs.size(); ++i) {
    double sum = 0.0;
    int count = 0;

    // Aggregate non-NaN distances
    for (size_t e = 0; e < embeddings.size(); ++e) {
      double distance = distances[i + e * libs.size()];
      if (!std::isnan(distance)) {
        sum += distance;
        ++count;
      }
    }

    // Compute mean, or set to NaN if all values are NaN
    if (count > 0) {
      row_means[i] = sum / count;
    } else {
      row_means[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  return row_means;
}

// Function to perform Simplex Projection Grid
double SimplexProjectionGrid(const std::vector<std::vector<std::vector<double>>>& embeddings,
                             const std::vector<double>& target, const std::vector<bool>& lib_indices,
                             const std::vector<bool>& pred_indices, int num_neighbors) {
  std::vector<double> pred(target.size(), std::numeric_limits<double>::quiet_NaN());

  for (size_t p = 0; p < pred_indices.size(); ++p) {
    if (!pred_indices[p]) continue;

    // Temporarily exclude the current prediction index from the library
    bool temp_lib = lib_indices[p];
    std::vector<bool> modified_lib_indices = lib_indices; // Create a copy of lib_indices
    modified_lib_indices[p] = false;

    // Get the indices of the library points
    std::vector<int> libs;
    for (size_t i = 0; i < modified_lib_indices.size(); ++i) {
      if (modified_lib_indices[i]) {
        libs.push_back(i);
      }
    }

    // Calculate distances
    std::vector<double> distances = DistCom(embeddings, libs, p);

    // Find the nearest neighbors, excluding NaN distances
    std::vector<int> neighbors;
    std::vector<double> valid_distances;
    for (size_t i = 0; i < distances.size(); ++i) {
      if (!std::isnan(distances[i])) {
        neighbors.push_back(i);
        valid_distances.push_back(distances[i]);
      }
    }

    // If fewer than num_neighbors valid distances, skip this prediction
    if (static_cast<int>(neighbors.size()) < num_neighbors) {
      continue; // Skip this prediction
    }

    // Sort valid distances and get the top num_neighbors
    std::vector<int> sorted_indices(neighbors.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&](int a, int b) {
      return valid_distances[a] < valid_distances[b];
    });

    // Get the top num_neighbors
    std::vector<int> top_neighbors;
    for (int i = 0; i < num_neighbors; ++i) {
      top_neighbors.push_back(neighbors[sorted_indices[i]]);
    }

    double min_distance = valid_distances[sorted_indices[0]];

    // Compute weights
    std::vector<double> weights(num_neighbors, 0.000001);
    if (min_distance == 0) {
      for (int i = 0; i < num_neighbors; ++i) {
        if (valid_distances[sorted_indices[i]] == 0) {
          weights[i] = 1.0;
        }
      }
    } else {
      for (int i = 0; i < num_neighbors; ++i) {
        weights[i] = std::exp(-valid_distances[sorted_indices[i]] / min_distance);
        if (weights[i] < 0.000001) {
          weights[i] = 0.000001;
        }
      }
    }

    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    // Make prediction
    double prediction = 0.0;
    for (int i = 0; i < num_neighbors; ++i) {
      prediction += weights[i] * target[libs[top_neighbors[i]]];
    }
    pred[p] = prediction / total_weight;
  }

  // Return the Pearson correlation between the target and the predicted values
  return PearsonCor(target, pred, true);
}

// GCCMSingle4Grid function
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b) {

  std::vector<std::pair<int, double>> x_xmap_y;

  for (int r = 1; r <= totalRow - lib_size + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size + 1; ++c) {

      // Initialize prediction and library indices
      std::vector<bool> pred_indices(totalRow * totalCol, false);
      std::vector<bool> lib_indices(totalRow * totalCol, false);

      // Set prediction indices
      for (const auto& p : pred) {
        pred_indices[locate(p.first, p.second, totalRow, totalCol)] = true;
      }

      // Exclude NA values in yPred from prediction indices
      for (size_t i = 0; i < yPred.size(); ++i) {
        if (std::isnan(yPred[i])) {
          pred_indices[i] = false;
        }
      }

      // Set library indices
      for (int i = r; i < r + lib_size; ++i) {
        for (int j = c; j < c + lib_size; ++j) {
          lib_indices[locate(i, j, totalRow, totalCol)] = true;
        }
      }

      // Check if more than half of the library is NA
      int na_count = 0;
      for (size_t i = 0; i < lib_indices.size(); ++i) {
        if (lib_indices[i] && std::isnan(yPred[i])) {
          ++na_count;
        }
      }
      if (na_count > (lib_size * lib_size) / 2) {
        continue;
      }

      // Run cross map and store results
      double results = SimplexProjectionGrid(xEmbedings, yPred, lib_indices, pred_indices, b);
      x_xmap_y.emplace_back(lib_size, results);
    }
  }

  return x_xmap_y;
}

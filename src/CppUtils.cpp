#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <limits>
#include "HelperFuns.h"

// Function to calculate the lagged indices
std::vector<std::vector<int>> CppLaggedIndices(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int lagNum) {
  int n = vec.size();
  std::vector<std::vector<int>> result(n);

  // Handle the case when lagNum is 0
  if (lagNum == 0) {
    for (int i = 0; i < n; ++i) {
      result[i] = {i};
    }
    return result;
  }

  // Handle the case when lagNum is greater than 0
  for (int i = 0; i < n; ++i) {
    std::unordered_set<int> visited;
    std::vector<int> current_neighbors;
    std::vector<int> next_neighbors;

    // Collect 1st level neighbors
    for (int j = 0; j < n; ++j) {
      if (nbmat[i][j] == 1 && i != j) {
        current_neighbors.push_back(j);
        visited.insert(j);
      }
    }

    // Collect neighbors up to lagNum
    for (int l = 1; l < lagNum; ++l) {
      for (int neighbor : current_neighbors) {
        for (int j = 0; j < n; ++j) {
          if (nbmat[neighbor][j] == 1 && i != j && visited.find(j) == visited.end()) {
            next_neighbors.push_back(j);
            visited.insert(j);
          }
        }
      }
      current_neighbors = next_neighbors;
      next_neighbors.clear();
    }

    // Convert set to vector and add to result
    result[i].insert(result[i].end(), visited.begin(), visited.end());

    // If no neighbors found, add NA
    if (result[i].empty()) {
      result[i].push_back(std::numeric_limits<int>::min());
    }
  }

  return result;
}

// Function to generate embeddings
std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int E) {
  int n = vec.size();
  std::vector<std::vector<double>> embeddings(n, std::vector<double>(E));

  for (int e = 0; e < E; ++e) {
    int lagNum = e;
    std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec, nbmat, lagNum);

    // Remove duplicates with previous lagNum
    if (e > 0) {
      std::vector<std::vector<int>> prev_lagged_indices = CppLaggedIndices(vec, nbmat, e - 1);
      for (int i = 0; i < n; ++i) {
        std::unordered_set<int> prev_set(prev_lagged_indices[i].begin(), prev_lagged_indices[i].end());
        std::vector<int> new_indices;
        for (int index : lagged_indices[i]) {
          if (prev_set.find(index) == prev_set.end()) {
            new_indices.push_back(index);
          }
        }
        lagged_indices[i] = new_indices;
        if (lagged_indices[i].empty()) {
          lagged_indices[i].push_back(std::numeric_limits<int>::min());
        }
      }
    }

    for (int i = 0; i < n; ++i) {
      std::vector<double> lagged_values;
      for (int index : lagged_indices[i]) {
        if (!checkIntNA(index)) {
          lagged_values.push_back(vec[index]);
        }
      }

      // Check if lagged_values is empty
      if (lagged_values.empty()) {
        embeddings[i][e] = std::numeric_limits<double>::quiet_NaN();
      } else {
        embeddings[i][e] = CppMean(lagged_values, true);
      }
      // embeddings[i][e] = CppMean(lagged_values, true);
    }
  }

  return embeddings;
}

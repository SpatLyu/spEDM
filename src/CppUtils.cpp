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

// // Function to calculate the pairwise absolute difference mean matrix
// std::vector<std::vector<double>> CppDist(const std::vector<std::vector<double>>& matrix) {
//   size_t n = matrix.size();
//   std::vector<std::vector<double>> result(n, std::vector<double>(n, std::numeric_limits<double>::quiet_NaN()));
//
//   for (size_t i = 0; i < n; ++i) {
//     for (size_t j = i; j < n; ++j) {
//       if (i == j) {
//         result[i][j] = 0.0; // Diagonal elements are zero
//       } else {
//         std::vector<double> abs_diff = CppAbs(matrix[i], matrix[j]);
//         result[i][j] = CppMean(abs_diff, true);
//         result[j][i] = result[i][j]; // Symmetric matrix
//       }
//     }
//   }
//
//   return result;
// }
//
// // Function to find the indices of the libsize + 1 closest elements for each row in distmat
// std::vector<std::vector<int>> CppClosestIndices(const std::vector<std::vector<double>>& distmat,
//                                                 int libsize) {
//   size_t n = distmat.size();
//   std::vector<std::vector<int>> closestIndices(n, std::vector<int>(libsize + 1));
//
//   for (size_t i = 0; i < n; ++i) {
//     // Create a vector to store the indices and distances
//     std::vector<std::pair<double, int>> distances;
//
//     for (size_t j = 0; j < n; ++j) {
//       if (i != j) { // Exclude the diagonal element
//         distances.push_back(std::make_pair(distmat[i][j], j));
//       }
//     }
//
//     // Sort the distances vector by the distance (first element of the pair)
//     std::sort(distances.begin(), distances.end(),
//               [](const std::pair<double, int>& a,
//                  const std::pair<double, int>& b) {
//       return a.first < b.first;
//     });
//
//     // Select the libsize + 1 closest indices
//     for (int k = 0; k < libsize + 1; ++k) {
//       closestIndices[i][k] = distances[k].second;
//     }
//   }
//
//   return closestIndices;
// }
//
// // Function to calculate the CCM weights
// std::vector<std::vector<double>> CCMWeight(const std::vector<std::vector<double>>& distmat,
//                                            const std::vector<std::vector<int>>& closestIndices,
//                                            int libsize) {
//   size_t n = distmat.size();
//   std::vector<std::vector<double>> result(n, std::vector<double>(libsize + 1));
//
//   for (size_t i = 0; i < n; ++i) {
//     std::vector<double> sdist;
//
//     // Extract the corresponding distances from distmat
//     for (int j = 0; j < libsize + 1; ++j) {
//       int index = closestIndices[i][j];
//       sdist.push_back(distmat[i][index]);
//     }
//
//     // Normalize sdist by dividing each element by the first element,
//     // taking the negative and applying the exponential function
//     double firstElement = sdist[0];
//     for (size_t j = 0; j < sdist.size(); ++j) {
//       sdist[j] = std::exp(-sdist[j] / firstElement);
//     }
//
//     // Normalize the resulting vector using CppSumNormalize with NA_rm set to true
//     std::vector<double> normalizedWeights = CppSumNormalize(sdist, true);
//
//     // Assign the normalized weights to the result
//     for (int j = 0; j < libsize + 1; ++j) {
//       result[i][j] = normalizedWeights[j];
//     }
//   }
//
//   return result;
// }
//
// // Function to perform Simplex Projection
// std::vector<double> SimplexProjection(const std::vector<double>& y,
//                                       const std::vector<double>& x,
//                                       const std::vector<std::vector<int>>& nbmat,
//                                       int libsize,
//                                       int E) {
//   // Generate embeddings for x
//   std::vector<std::vector<double>> singleE = GenEmbeddings(x, nbmat, E);
//
//   // Calculate the distance matrix for singleE
//   std::vector<std::vector<double>> singleDist = CppDist(singleE);
//
//   // Find the closest indices for each row in singleDist
//   std::vector<std::vector<int>> singleCI = CppClosestIndices(singleDist, libsize);
//
//   // Calculate the weights for each row in singleCI
//   std::vector<std::vector<double>> singleW = CCMWeight(singleDist, singleCI, libsize);
//
//   // Construct singleY by extracting values from y based on singleCI
//   size_t n = y.size();
//   std::vector<std::vector<double>> singleY(n, std::vector<double>(libsize + 1));
//   for (size_t i = 0; i < n; ++i) {
//     for (int j = 0; j < libsize + 1; ++j) {
//       int index = singleCI[i][j];
//       singleY[i][j] = y[index];
//     }
//   }
//
//   // Calculate the result by multiplying singleY and singleW and summing the results
//   std::vector<double> result(n);
//   for (size_t i = 0; i < n; ++i) {
//     std::vector<double> sumvec;
//     for (int j = 0; j < libsize + 1; ++j) {
//       sumvec.push_back(singleY[i][j] * singleW[i][j]);
//     }
//     result[i] = CppSum(sumvec, true);
//   }
//
//   return result;
// }

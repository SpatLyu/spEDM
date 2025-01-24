#include <iostream>
#include <vector>
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <limits> // for std::numeric_limits

// Function to compute lagged neighborhoods for a given lag number
// **Note that the return value corresponds to the cumulative neighbor indices for lagNum**
std::vector<std::vector<int>> CppLaggedVar4Lattice(std::vector<std::vector<int>> spNeighbor,
                                                   int lagNum) {
  // If lagNum is less than 1, return an empty vector
  if (lagNum < 1) {
    return {};
  }

  // Initialize the lagged neighborhood as a copy of spNeighbor
  std::vector<std::vector<int>> lagSpNeighbor = spNeighbor;

  // If lagNum is greater than 1, expand the neighborhoods
  if (lagNum > 1) {
    std::vector<std::vector<int>> curSpNeighbor = spNeighbor;

    // Iterate from 1 to lagNum - 1 to expand the neighborhoods
    for (int lag = 1; lag < lagNum; ++lag) {
      std::vector<std::vector<int>> preSpNeighbor = curSpNeighbor;

      // Update the current neighborhood for each node
      for (size_t i = 0; i < preSpNeighbor.size(); ++i) {
        std::vector<int> curChain = preSpNeighbor[i];
        std::vector<int> newRings = curChain;

        // Expand the neighborhood by including neighbors of neighbors
        for (int neigh : curChain) {
          if (neigh >= 0) {
            std::vector<int> nextChain = spNeighbor[neigh]; // Use 0-based index
            newRings.insert(newRings.end(), nextChain.begin(), nextChain.end());
          }
        }

        // Remove duplicates and sort the new neighborhood
        std::sort(newRings.begin(), newRings.end());
        newRings.erase(std::unique(newRings.begin(), newRings.end()), newRings.end());

        // Update the current neighborhood
        curSpNeighbor[i] = newRings;
      }
    }

    // Remove the original neighbors and the node itself from the lagged neighborhood
    for (size_t i = 0; i < curSpNeighbor.size(); ++i) {
      std::vector<int> newRings = curSpNeighbor[i];
      std::vector<int> original = spNeighbor[i];
      original.push_back(i); // Add the node itself (0-based index)

      // Remove original neighbors and the node itself
      std::vector<int> filteredRings;
      for (int ring : newRings) {
        if (std::find(original.begin(), original.end(), ring) == original.end()) {
          filteredRings.push_back(ring);
        }
      }

      // Update the lagged neighborhood
      lagSpNeighbor[i] = filteredRings;
    }
  }

  return lagSpNeighbor;
}

// Function to generate embeddings for a given vector and neighborhood matrix
std::vector<std::vector<double>> GenLatticeEmbeddings(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int E) {
  // Get the number of nodes
  int n = vec.size();

  // Initialize the embeddings matrix with NaN values
  std::vector<std::vector<double>> xEmbedings(n, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));

  // Compute embeddings for each lag number from 1 to E
  for (int lagNum = 1; lagNum <= E; ++lagNum) {
    // Compute the lagged neighborhoods
    std::vector<std::vector<int>> laggedResults = CppLaggedVar4Lattice(nb, lagNum);

    // Remove duplicates with previous lagNum (if lagNum >= 2)
    if (lagNum >= 2) {
      std::vector<std::vector<int>> prev_laggedResults = CppLaggedVar4Lattice(nb, lagNum - 1);
      for (int i = 0; i < n; ++i) {
        // Convert previous lagged results to a set for fast lookup
        std::unordered_set<int> prev_set(prev_laggedResults[i].begin(), prev_laggedResults[i].end());

        // Remove duplicates from current lagged results
        std::vector<int> new_indices;
        for (int index : laggedResults[i]) {
          if (prev_set.find(index) == prev_set.end()) {
            new_indices.push_back(index);
          }
        }

        // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
        if (new_indices.empty()) {
          new_indices.push_back(std::numeric_limits<int>::min());
        }

        // Update the lagged results
        laggedResults[i] = new_indices;
      }
    }

    // Compute the mean of neighbor values for each node
    for (size_t l = 0; l < laggedResults.size(); ++l) {
      std::vector<int> neighbors = laggedResults[l];

      // If the neighbors are empty or contain only the special value, leave the embedding as NaN
      if (neighbors.empty() || (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min())) {
        continue;
      }

      // Compute the mean of neighbor values
      double sum = std::accumulate(neighbors.begin(), neighbors.end(), 0.0, [&](double acc, int idx) {
        return acc + vec[idx];
      });
      xEmbedings[l][lagNum - 1] = sum / neighbors.size();
    }
  }

  return xEmbedings;
}

// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include <unordered_set>
// #include <limits>
//
// // Function to calculate the lagged indices
// std::vector<std::vector<int>> CppLaggedIndices(const std::vector<double>& vec,
//                                                const std::vector<std::vector<int>>& nbmat,
//                                                int lagNum) {
//   int n = vec.size();
//   std::vector<std::vector<int>> result(n);
//
//   // Handle the case when lagNum is 0
//   if (lagNum == 0) {
//     for (int i = 0; i < n; ++i) {
//       result[i] = {i};
//     }
//     return result;
//   }
//
//   // Handle the case when lagNum is greater than 0
//   for (int i = 0; i < n; ++i) {
//     std::unordered_set<int> visited;
//     std::vector<int> current_neighbors;
//     std::vector<int> next_neighbors;
//
//     // Collect 1st level neighbors
//     for (int j = 0; j < n; ++j) {
//       if (nbmat[i][j] == 1 && i != j) {
//         current_neighbors.push_back(j);
//         visited.insert(j);
//       }
//     }
//
//     // Collect neighbors up to lagNum
//     for (int l = 1; l < lagNum; ++l) {
//       for (int neighbor : current_neighbors) {
//         for (int j = 0; j < n; ++j) {
//           if (nbmat[neighbor][j] == 1 && i != j && visited.find(j) == visited.end()) {
//             next_neighbors.push_back(j);
//             visited.insert(j);
//           }
//         }
//       }
//       current_neighbors = next_neighbors;
//       next_neighbors.clear();
//     }
//
//     // Convert set to vector and add to result
//     result[i].insert(result[i].end(), visited.begin(), visited.end());
//
//     // If no neighbors found, add NA
//     if (result[i].empty()) {
//       result[i].push_back(std::numeric_limits<int>::min());
//     }
//   }
//
//   return result;
// }
//
// // Function to generate embeddings
// std::vector<std::vector<double>> GenLatticeEmbeddings(const std::vector<double>& vec,
//                                                       const std::vector<std::vector<int>>& nbmat,
//                                                       int E) {
//   int n = vec.size();
//   std::vector<std::vector<double>> embeddings(n, std::vector<double>(E));
//
//   for (int e = 0; e < E; ++e) {
//     int lagNum = e;
//     std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec, nbmat, lagNum);
//
//     // Remove duplicates with previous lagNum
//     if (e > 0) {
//       std::vector<std::vector<int>> prev_lagged_indices = CppLaggedIndices(vec, nbmat, e - 1);
//       for (int i = 0; i < n; ++i) {
//         std::unordered_set<int> prev_set(prev_lagged_indices[i].begin(), prev_lagged_indices[i].end());
//         std::vector<int> new_indices;
//         for (int index : lagged_indices[i]) {
//           if (prev_set.find(index) == prev_set.end()) {
//             new_indices.push_back(index);
//           }
//         }
//         lagged_indices[i] = new_indices;
//         if (lagged_indices[i].empty()) {
//           lagged_indices[i].push_back(std::numeric_limits<int>::min());
//         }
//       }
//     }
//
//     for (int i = 0; i < n; ++i) {
//       std::vector<double> lagged_values;
//       for (int index : lagged_indices[i]) {
//         if (!checkIntNA(index)) {
//           lagged_values.push_back(vec[index]);
//         }
//       }
//
//       // Check if lagged_values is empty
//       if (lagged_values.empty()) {
//         embeddings[i][e] = std::numeric_limits<double>::quiet_NaN();
//       } else {
//         embeddings[i][e] = CppMean(lagged_values, true);
//       }
//       // embeddings[i][e] = CppMean(lagged_values, true);
//     }
//   }
//
//   return embeddings;
// }

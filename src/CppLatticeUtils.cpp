#include <iostream>
#include <vector>
#include <algorithm> // for std::sort and std::unique
#include <numeric>   // for std::accumulate
#include <Rcpp.h>

std::vector<std::vector<int>> nb2vec(Rcpp::List nb) {
  // Get the number of elements in the nb object
  int n = nb.size();

  // Create a vector<vector<int>> to store the result
  std::vector<std::vector<int>> result(n);

  // Iterate over each element in the nb object
  for (int i = 0; i < n; ++i) {
    // Get the current element (should be an integer vector)
    Rcpp::IntegerVector current_nb = nb[i];

    // Create a vector<int> to store the current subset of elements
    std::vector<int> current_subset;

    // Iterate over each element in the current subset
    for (int j = 0; j < current_nb.size(); ++j) {
      // Subtract one from each element to convert from R's 1-based indexing to C++'s 0-based indexing
      current_subset.push_back(current_nb[j] - 1);
    }

    // Add the current subset to the result
    result[i] = current_subset;
  }

  return result;
}

// Function to compute lagged neighborhoods for a given lag number
std::vector<std::vector<int>> CppLaggedVar4Lattice(std::vector<std::vector<int>> spNeighbor,
                                                   int lagNum) {
  // If lagNum is less than 1, return an empty vector
  if (lagNum < 1) {
    return {};
  }

  // If lagNum is 1, return a copy of spNeighbor
  if (lagNum == 1) {
    return spNeighbor;
  }

  // Initialize the current neighborhood as a copy of spNeighbor
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
        if (neigh > 0) {
          std::vector<int> nextChain = spNeighbor[neigh - 1]; // Convert to 0-based index
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
  std::vector<std::vector<int>> lagSpNeighbor = curSpNeighbor;
  for (size_t i = 0; i < curSpNeighbor.size(); ++i) {
    std::vector<int> newRings = curSpNeighbor[i];
    std::vector<int> original = spNeighbor[i];
    original.push_back(i + 1); // Add the node itself (convert to 1-based index)

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

  return lagSpNeighbor;
}

// Function to generate embeddings for a given vector and neighborhood matrix
std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
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

    // Compute the mean of neighbor values for each node
    for (size_t l = 0; l < laggedResults.size(); ++l) {
      std::vector<int> neighbors = laggedResults[l];

      // Convert neighbors to 0-based index
      for (int& neighbor : neighbors) {
        neighbor -= 1;
      }

      // Compute the mean of neighbor values
      if (!neighbors.empty()) {
        double sum = std::accumulate(neighbors.begin(), neighbors.end(), 0.0, [&](double acc, int idx) {
          return acc + vec[idx];
        });
        xEmbedings[l][lagNum - 1] = sum / neighbors.size();
      }
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
// std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
//                                                const std::vector<std::vector<int>>& nbmat,
//                                                int E) {
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

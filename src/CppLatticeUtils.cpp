#include <iostream>
#include <vector>
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <limits> // for std::numeric_limits
#include <cmath> // For std::isnan

/**
 * Computes lagged neighborhoods for a given lag number, expanding the neighborhoods iteratively
 * by including neighbors of neighbors up to the specified lag number.
 *
 * Parameters:
 *   spNeighbor - A 2D vector representing the spatial neighbors for each spatial unit, where each element is a list of neighbors.
 *   lagNum     - The number of lags to expand the neighborhoods.
 *                A lagNum of 0 means the original spNeighbor is returned.
 *                A lagNum of 1 means only the immediate neighbors are considered.
 *                A lagNum >= 2 means the neighborhoods are expanded by including neighbors of neighbors, and the result is filtered to remove the immediate neighbors.
 *
 * Returns:
 *   A 2D vector where each element is a list of cumulative neighbor indices for a given spatial unit,
 *   including neighbors up to the specified lagNum. If lagNum is 0, the original spNeighbor is returned.
 *   If lagNum >= 2, the result is filtered to remove the immediate neighbors (lagNum = 1). If the filtered result is empty, it is filled with std::numeric_limits<int>::min().
 */
std::vector<std::vector<int>> CppLaggedNeighbor4Lattice(std::vector<std::vector<int>> spNeighbor,
                                                        int lagNum) {
  // If lagNum is less than 0, return an empty vector
  if (lagNum < 0) {
    return {};
  }

  // If lagNum is 0, return the original spNeighbor
  if (lagNum == 0) {
    return spNeighbor;
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

      // If the filtered result is empty, fill it with std::numeric_limits<int>::min()
      if (filteredRings.empty()) {
        filteredRings.push_back(std::numeric_limits<int>::min());
      }

      // Update the lagged neighborhood
      lagSpNeighbor[i] = filteredRings;
    }
  }

  return lagSpNeighbor;
}

/**
 * Computes lagged neighborhoods for a given lag number, expanding the neighborhoods iteratively
 * by including neighbors of neighbors up to the specified lag number.
 *
 * Parameters:
 *   spNeighbor - A 2D vector representing the spatial neighbors for each spatial unit, where each element is a list of neighbors.
 *   lagNum     - The number of lags to expand the neighborhoods.
 *                A lagNum of 1 means only the immediate neighbors are considered.
 *
 * Returns:
 *   A 2D vector where each element is a list of cumulative neighbor indices for a given spatial unit,
 *   including neighbors up to the specified lagNum. If lagNum is less than 1, an empty vector is returned.
 *
 * Note:
 *   The return value corresponds to the cumulative neighbor indices for the specified lagNum.
 *   The neighborhoods are expanded by including neighbors of neighbors, and duplicates are removed at each step.
 */
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

/**
 * Generates embeddings for a given vector and neighborhood matrix by computing the mean of neighbor values
 * for each spatial unit, considering both the immediate neighbors and neighbors up to a specified lag number.
 *
 * Parameters:
 *   vec  - A vector of values, one for each spatial unit, to be embedded.
 *   nb   - A 2D matrix where each row represents the neighbors of a spatial unit.
 *   E    - The embedding dimension, specifying how many lags to consider in the embeddings.
 *   tau  - The spatial lag step for constructing lagged state-space vectors.
 *
 * Returns:
 *   A 2D vector representing the embeddings for each spatial unit, where each spatial unit has a row in the matrix
 *   corresponding to its embedding values for each lag number. If no valid embedding columns remain after removing
 *   columns containing only NaN values, a filtered matrix is returned.
 */
std::vector<std::vector<double>> GenLatticeEmbeddings(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int E,
    int tau)
{
  // Get the number of nodes
  int n = vec.size();

  // Initialize the embeddings matrix with NaN values
  std::vector<std::vector<double>> xEmbedings(n, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));

  if (tau == 0) {
    // When E >= 1, fill the first column of xEmbedings with the values from vec
    if (E >= 1) {
      for (int i = 0; i < n; ++i) {
        xEmbedings[i][0] = vec[i]; // Fill the first column with vec values
      }
    }

    // Compute embeddings for each lag number from 1 to E
    for (int lagNum = 1; lagNum < E; ++lagNum) {
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
        xEmbedings[l][lagNum] = sum / neighbors.size();
      }
    }
  } else {
    // Compute embeddings for each lag number from tau to E+tau-1
    for (int lagNum = tau; lagNum < E + tau; ++lagNum) {
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
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns; // To store indices of valid columns

  // Iterate over each column to check if it contains any non-NaN values
  for (size_t col = 0; col < xEmbedings[0].size(); ++col) {
    bool isAllNaN = true;
    for (size_t row = 0; row < xEmbedings.size(); ++row) {
      if (!std::isnan(xEmbedings[row][col])) {
        isAllNaN = false;
        break;
      }
    }
    if (!isAllNaN) {
      validColumns.push_back(col); // Store the index of valid columns
    }
  }

  // If no columns are removed, return the original xEmbedings
  if (validColumns.size() == xEmbedings[0].size()) {
    return xEmbedings;
  } else {
    // // Issue a warning if any columns are removed
    // std::cerr << "Warning: remove all-NA embedding vector columns caused by excessive embedding dimension E selection." << std::endl;

    // Construct the filtered embeddings matrix
    std::vector<std::vector<double>> filteredEmbeddings;
    for (size_t row = 0; row < xEmbedings.size(); ++row) {
      std::vector<double> filteredRow;
      for (size_t col : validColumns) {
        filteredRow.push_back(xEmbedings[row][col]);
      }
      filteredEmbeddings.push_back(filteredRow);
    }

    // Return the filtered embeddings matrix
    return filteredEmbeddings;
  }
}

// #include <iostream>
// #include <vector>
// #include <numeric>   // for std::accumulate
// #include <algorithm> // for std::sort, std::unique, std::accumulate
// #include <unordered_set> // for std::unordered_set
// #include <limits> // for std::numeric_limits
// #include <cmath> // For std::isnan
//
// /**
//  * Computes lagged neighborhoods for a given lag number, expanding the neighborhoods iteratively
//  * by including neighbors of neighbors up to the specified lag number.
//  *
//  * Parameters:
//  *   spNeighbor - A 2D vector representing the spatial neighbors for each spatial unit, where each element is a list of neighbors.
//  *   lagNum     - The number of lags to expand the neighborhoods.
//  *                A lagNum of 1 means only the immediate neighbors are considered.
//  *
//  * Returns:
//  *   A 2D vector where each element is a list of cumulative neighbor indices for a given spatial unit,
//  *   including neighbors up to the specified lagNum. If lagNum is less than 1, an empty vector is returned.
//  *
//  * Note:
//  *   The return value corresponds to the cumulative neighbor indices for the specified lagNum.
//  *   The neighborhoods are expanded by including neighbors of neighbors, and duplicates are removed at each step.
//  */
// std::vector<std::vector<int>> CppLaggedVar4Lattice(std::vector<std::vector<int>> spNeighbor,
//                                                    int lagNum) {
//   // If lagNum is less than 1, return an empty vector
//   if (lagNum < 1) {
//     return {};
//   }
//
//   // Initialize the lagged neighborhood as a copy of spNeighbor
//   std::vector<std::vector<int>> lagSpNeighbor = spNeighbor;
//
//   // If lagNum is greater than 1, expand the neighborhoods
//   if (lagNum > 1) {
//     std::vector<std::vector<int>> curSpNeighbor = spNeighbor;
//
//     // Iterate from 1 to lagNum - 1 to expand the neighborhoods
//     for (int lag = 1; lag < lagNum; ++lag) {
//       std::vector<std::vector<int>> preSpNeighbor = curSpNeighbor;
//
//       // Update the current neighborhood for each node
//       for (size_t i = 0; i < preSpNeighbor.size(); ++i) {
//         std::vector<int> curChain = preSpNeighbor[i];
//         std::vector<int> newRings = curChain;
//
//         // Expand the neighborhood by including neighbors of neighbors
//         for (int neigh : curChain) {
//           if (neigh >= 0) {
//             std::vector<int> nextChain = spNeighbor[neigh]; // Use 0-based index
//             newRings.insert(newRings.end(), nextChain.begin(), nextChain.end());
//           }
//         }
//
//         // Remove duplicates and sort the new neighborhood
//         std::sort(newRings.begin(), newRings.end());
//         newRings.erase(std::unique(newRings.begin(), newRings.end()), newRings.end());
//
//         // Update the current neighborhood
//         curSpNeighbor[i] = newRings;
//       }
//     }
//
//     // Remove the original neighbors and the node itself from the lagged neighborhood
//     for (size_t i = 0; i < curSpNeighbor.size(); ++i) {
//       std::vector<int> newRings = curSpNeighbor[i];
//       std::vector<int> original = spNeighbor[i];
//       original.push_back(i); // Add the node itself (0-based index)
//
//       // Remove original neighbors and the node itself
//       std::vector<int> filteredRings;
//       for (int ring : newRings) {
//         if (std::find(original.begin(), original.end(), ring) == original.end()) {
//           filteredRings.push_back(ring);
//         }
//       }
//
//       // Update the lagged neighborhood
//       lagSpNeighbor[i] = filteredRings;
//     }
//   }
//
//   return lagSpNeighbor;
// }
//
// /**
//  * Generates embeddings for a given vector and neighborhood matrix by computing the mean of neighbor values
//  * for each spatial unit, considering both the immediate neighbors and neighbors up to a specified lag number.
//  *
//  * Parameters:
//  *   vec  - A vector of values, one for each spatial unit, to be embedded.
//  *   nb   - A 2D matrix where each row represents the neighbors of a spatial unit.
//  *   E    - The embedding dimension, specifying how many lags to consider in the embeddings.
//  *   tau  - The spatial lag step for constructing lagged state-space vectors.
//  *
//  * Returns:
//  *   A 2D vector representing the embeddings for each spatial unit, where each spatial unit has a row in the matrix
//  *   corresponding to its embedding values for each lag number. If no valid embedding columns remain after removing
//  *   columns containing only NaN values, a filtered matrix is returned.
//  */
// std::vector<std::vector<double>> GenLatticeEmbeddings(
//     const std::vector<double>& vec,
//     const std::vector<std::vector<int>>& nb,
//     int E,
//     int tau)
// {
//   // Get the number of nodes
//   int n = vec.size();
//
//   // Initialize the embeddings matrix with NaN values
//   std::vector<std::vector<double>> xEmbedings(n, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));
//
//   if (tau == 0) {
//     // When E >= 1, fill the first column of xEmbedings with the values from vec
//     if (E >= 1) {
//       for (int i = 0; i < n; ++i) {
//         xEmbedings[i][0] = vec[i]; // Fill the first column with vec values
//       }
//     }
//
//     // Compute embeddings for each lag number from 1 to E
//     for (int lagNum = 1; lagNum < E; ++lagNum) {
//       // Compute the lagged neighborhoods
//       std::vector<std::vector<int>> laggedResults = CppLaggedVar4Lattice(nb, lagNum);
//
//       // Remove duplicates with previous lagNum (if lagNum >= 2)
//       if (lagNum >= 2) {
//         std::vector<std::vector<int>> prev_laggedResults = CppLaggedVar4Lattice(nb, lagNum - 1);
//         for (int i = 0; i < n; ++i) {
//           // Convert previous lagged results to a set for fast lookup
//           std::unordered_set<int> prev_set(prev_laggedResults[i].begin(), prev_laggedResults[i].end());
//
//           // Remove duplicates from current lagged results
//           std::vector<int> new_indices;
//           for (int index : laggedResults[i]) {
//             if (prev_set.find(index) == prev_set.end()) {
//               new_indices.push_back(index);
//             }
//           }
//
//           // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
//           if (new_indices.empty()) {
//             new_indices.push_back(std::numeric_limits<int>::min());
//           }
//
//           // Update the lagged results
//           laggedResults[i] = new_indices;
//         }
//       }
//
//       // Compute the mean of neighbor values for each node
//       for (size_t l = 0; l < laggedResults.size(); ++l) {
//         std::vector<int> neighbors = laggedResults[l];
//
//         // If the neighbors are empty or contain only the special value, leave the embedding as NaN
//         if (neighbors.empty() || (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min())) {
//           continue;
//         }
//
//         // Compute the mean of neighbor values
//         double sum = std::accumulate(neighbors.begin(), neighbors.end(), 0.0, [&](double acc, int idx) {
//           return acc + vec[idx];
//         });
//         xEmbedings[l][lagNum] = sum / neighbors.size();
//       }
//     }
//   } else {
//     // Compute embeddings for each lag number from tau to E+tau-1
//     for (int lagNum = tau; lagNum < E + tau; ++lagNum) {
//       // Compute the lagged neighborhoods
//       std::vector<std::vector<int>> laggedResults = CppLaggedVar4Lattice(nb, lagNum);
//
//       // Remove duplicates with previous lagNum (if lagNum >= 2)
//       if (lagNum >= 2) {
//         std::vector<std::vector<int>> prev_laggedResults = CppLaggedVar4Lattice(nb, lagNum - 1);
//         for (int i = 0; i < n; ++i) {
//           // Convert previous lagged results to a set for fast lookup
//           std::unordered_set<int> prev_set(prev_laggedResults[i].begin(), prev_laggedResults[i].end());
//
//           // Remove duplicates from current lagged results
//           std::vector<int> new_indices;
//           for (int index : laggedResults[i]) {
//             if (prev_set.find(index) == prev_set.end()) {
//               new_indices.push_back(index);
//             }
//           }
//
//           // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
//           if (new_indices.empty()) {
//             new_indices.push_back(std::numeric_limits<int>::min());
//           }
//
//           // Update the lagged results
//           laggedResults[i] = new_indices;
//         }
//       }
//
//       // Compute the mean of neighbor values for each node
//       for (size_t l = 0; l < laggedResults.size(); ++l) {
//         std::vector<int> neighbors = laggedResults[l];
//
//         // If the neighbors are empty or contain only the special value, leave the embedding as NaN
//         if (neighbors.empty() || (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min())) {
//           continue;
//         }
//
//         // Compute the mean of neighbor values
//         double sum = std::accumulate(neighbors.begin(), neighbors.end(), 0.0, [&](double acc, int idx) {
//           return acc + vec[idx];
//         });
//         xEmbedings[l][lagNum - 1] = sum / neighbors.size();
//       }
//     }
//   }
//
//   // Calculate validColumns (indices of columns that are not entirely NaN)
//   std::vector<size_t> validColumns; // To store indices of valid columns
//
//   // Iterate over each column to check if it contains any non-NaN values
//   for (size_t col = 0; col < xEmbedings[0].size(); ++col) {
//     bool isAllNaN = true;
//     for (size_t row = 0; row < xEmbedings.size(); ++row) {
//       if (!std::isnan(xEmbedings[row][col])) {
//         isAllNaN = false;
//         break;
//       }
//     }
//     if (!isAllNaN) {
//       validColumns.push_back(col); // Store the index of valid columns
//     }
//   }
//
//   // If no columns are removed, return the original xEmbedings
//   if (validColumns.size() == xEmbedings[0].size()) {
//     return xEmbedings;
//   } else {
//     // // Issue a warning if any columns are removed
//     // std::cerr << "Warning: remove all-NA embedding vector columns caused by excessive embedding dimension E selection." << std::endl;
//
//     // Construct the filtered embeddings matrix
//     std::vector<std::vector<double>> filteredEmbeddings;
//     for (size_t row = 0; row < xEmbedings.size(); ++row) {
//       std::vector<double> filteredRow;
//       for (size_t col : validColumns) {
//         filteredRow.push_back(xEmbedings[row][col]);
//       }
//       filteredEmbeddings.push_back(filteredRow);
//     }
//
//     // Return the filtered embeddings matrix
//     return filteredEmbeddings;
//   }
// }

#include <iostream>
#include <vector>
#include <queue> // for std::queue
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <unordered_map> // for std::unordered_map
#include <limits> // for std::numeric_limits
#include <cmath> // For std::isnan

/**
 * Computes the lagged neighbors for a lattice structure up to a specified lag number.
 * This function recursively expands the neighbors at each lag step, starting with direct neighbors
 * (lag 0), and including neighbors from previous lags, until reaching the specified lag number.
 *
 * For lagNum = 0, each spatial unit is its own neighbor.
 * For lagNum >= 1, the function accumulates neighbors from all previous lags and deduplicates the results.
 * Empty results are filled with `std::numeric_limits<int>::min()` to indicate no neighbors.
 *
 * Parameters:
 *   spNeighbor - A 2D vector where each element contains indices of immediate neighbors for each spatial unit.
 *   lagNum     - The number of lag steps to compute (must be non-negative).
 *
 * Returns:
 *   A 2D vector where each element represents the list of lagged neighbors for a spatial unit.
 */
std::vector<std::vector<int>> CppLaggedNeighbor4Lattice(const std::vector<std::vector<int>>& spNeighbor,
                                                        int lagNum) {
  // Handle negative lagNum: return empty vector
  if (lagNum < 0) {
    return {};
  }

  // If lagNum is 0, return a vector of indices
  if (lagNum == 0) {
    std::vector<std::vector<int>> result;
    for (size_t i = 0; i < spNeighbor.size(); ++i) {
      result.push_back({static_cast<int>(i)});
    }
    return result;
  }

  // // Handle lagNum=1: return the immediate neighbors directly
  // if (lagNum == 1) {
  //   return spNeighbor;
  // }

  // Recursively compute results for lagNum-1
  std::vector<std::vector<int>> prevResult = CppLaggedNeighbor4Lattice(spNeighbor, lagNum - 1);
  std::vector<std::vector<int>> currentResult;

  int n = spNeighbor.size();
  // Process each spatial unit to compute current lagNum's neighbors
  for (int i = 0; i < n; ++i) {
    // Check if prevResult[i] size is equal to n
    if (prevResult[i].size() == spNeighbor.size()) {
      currentResult.push_back(prevResult[i]);
      continue; // Skip further processing for this index
    }

    std::unordered_set<int> mergedSet;

    // Add previous lag results (lag from 0 to lagNum-1)
    for (int elem : prevResult[i]) {
      if (elem != std::numeric_limits<int>::min()) {
        mergedSet.insert(elem);
      }
    }

    // Collect new elements from neighbors of previous lag's results
    std::unordered_set<int> newElements;
    for (int j : prevResult[i]) {
      // Skip invalid indices and placeholder min value
      if (j == std::numeric_limits<int>::min() || j < 0 || j >= n) {
        continue;
      }
      // Aggregate neighbors of j from spNeighbor
      for (int k : spNeighbor[j]) {
        newElements.insert(k);
      }
    }

    // Merge new elements into the set
    for (int elem : newElements) {
      mergedSet.insert(elem);
    }

    // Convert set to sorted vector and deduplicate
    std::vector<int> vec(mergedSet.begin(), mergedSet.end());
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    // Handle empty result by filling with min value
    if (vec.empty()) {
      vec.push_back(std::numeric_limits<int>::min());
    }

    currentResult.push_back(vec);
  }

  return currentResult;
}

/**
 * Computes the lagged values for a given vector based on the neighborhood structure and lag number.
 * This function first determines the lagged neighbors for each spatial unit using
 * the `CppLaggedNeighbor4Lattice` function. If `lagNum > 0`, it removes duplicate indices that
 * appeared in previous lag levels to ensure each lag level captures only new neighbors.
 *
 * For each spatial unit, the function extracts values from `vec` corresponding to the computed
 * lagged neighbors. If no valid neighbors exist, the function fills the result with `NaN`.
 *
 * Parameters:
 *   vec    - A vector of double values representing the spatial data for each unit.
 *   nb     - A 2D vector where each row contains indices of immediate neighbors in the lattice.
 *   lagNum - The number of lag steps to compute (must be non-negative).
 *
 * Returns:
 *   A 2D vector where each element contains the lagged values corresponding to the computed
 *   lagged neighbors for each spatial unit.
 */
std::vector<std::vector<double>> CppLaggedVar4Lattice(const std::vector<double>& vec,
                                                      const std::vector<std::vector<int>>& nb,
                                                      int lagNum) {
  int n = vec.size();

  // Compute the lagged neighbors using the provided function
  std::vector<std::vector<int>> laggedNeighbors = CppLaggedNeighbor4Lattice(nb, lagNum);
  // Remove duplicates with previous lagNum (if lagNum > 0)
  if (lagNum > 0) {
    std::vector<std::vector<int>> prevLaggedResults = CppLaggedNeighbor4Lattice(nb, lagNum - 1);
    for (int i = 0; i < n; ++i) {
      // Convert previous lagged results to a set for fast lookup
      std::unordered_set<int> prevSet(prevLaggedResults[i].begin(), prevLaggedResults[i].end());

      // Remove duplicates from current lagged results
      std::vector<int> newIndices;
      for (int index : laggedNeighbors[i]) {
        if (prevSet.find(index) == prevSet.end()) {
          newIndices.push_back(index);
        }
      }

      // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
      if (newIndices.empty()) {
        newIndices.push_back(std::numeric_limits<int>::min());
      }

      // Update the lagged results
      laggedNeighbors[i] = newIndices;
    }
  }

  // Initialize the result vector with the same number of rows as the lagged neighbors
  std::vector<std::vector<double>> result(laggedNeighbors.size());

  // Iterate over each point in the lattice
  for (size_t i = 0; i < laggedNeighbors.size(); ++i) {
    // Initialize the lagged values for the current point
    std::vector<double> laggedValues;

    if (laggedNeighbors[i].size() == 1 && laggedNeighbors[i][0] == std::numeric_limits<int>::min()){
      // If the index is out of bounds, push a default value (e.g., nan)
      laggedValues.push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      // Iterate over each neighbor index and extract the corresponding value from `vec`
      for (int neighborIndex : laggedNeighbors[i]) {
        // Check if the neighbor index is valid
        if (neighborIndex >= 0 && neighborIndex < n) {
          laggedValues.push_back(vec[neighborIndex]);
        } else {
          // If the index is out of bounds, push a default value (e.g., nan)
          laggedValues.push_back(std::numeric_limits<double>::quiet_NaN());
        }
      }
    }

    // Add the lagged values to the result
    result[i] = laggedValues;
  }

  return result;
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

  // Precompute lagged neighbor results for all required lagNum values
  std::unordered_map<int, std::vector<std::vector<int>>> laggedResultsMap;

  // Determine the range of lagNum values based on tau
  int startLagNum = (tau == 0) ? 0 : tau;
  int endLagNum = (tau == 0) ? E - 1 : E * tau;
  int step = (tau == 0) ? 1 : tau;

  for (int lagNum = 0; lagNum <= endLagNum; ++lagNum){
    if (lagNum == 0) { // return the current index (C++ based 0 index) for each spatial unit;
      std::vector<std::vector<int>> result_temp;
      for (size_t i = 0; i < nb.size(); ++i) {
        result_temp.push_back({static_cast<int>(i)});
      }
      laggedResultsMap[lagNum] = result_temp;
    } else { // when lagNum > 0, recursively compute results for lagNum-1;
      std::vector<std::vector<int>> prevResult = laggedResultsMap[lagNum - 1];
      std::vector<std::vector<int>> currentResult;

      // Process each spatial unit to compute current lagNum's neighbors
      for (int i = 0; i < n; ++i) {
        // Check if prevResult[i] size is equal to n
        if (prevResult[i].size() == nb.size()) {
          currentResult.push_back(prevResult[i]);
          continue; // Skip further processing for this index
        }

        std::unordered_set<int> mergedSet;

        // Add previous lag results (lag from 0 to lagNum-1)
        for (int elem : prevResult[i]) {
          if (elem != std::numeric_limits<int>::min()) {
            mergedSet.insert(elem);
          }
        }

        // Collect new elements from neighbors of previous lag's results
        std::unordered_set<int> newElements;
        for (int j : prevResult[i]) {
          // Skip invalid indices and placeholder min value
          if (j == std::numeric_limits<int>::min() || j < 0 || j >= n) {
            continue;
          }
          // Aggregate neighbors of j from spNeighbor
          for (int k : nb[j]) {
            newElements.insert(k);
          }
        }

        // Merge new elements into the set
        for (int elem : newElements) {
          mergedSet.insert(elem);
        }

        // Convert set to sorted vector and deduplicate
        std::vector<int> vec(mergedSet.begin(), mergedSet.end());
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

        // Handle empty result by filling with min value
        if (vec.empty()) {
          vec.push_back(std::numeric_limits<int>::min());
        }

        currentResult.push_back(vec);
      }

      laggedResultsMap[lagNum] = currentResult;
    }
  }

  // // Generate a sequence of lagNum values that need to be computed, including lagNum and lagNum - step
  // std::unordered_set<int> lagNumNeedSet;
  //
  // for (int lagNum = startLagNum; lagNum <= endLagNum; lagNum += step) {
  //   lagNumNeedSet.insert(lagNum);
  //   if (lagNum > 0) {
  //     lagNumNeedSet.insert(lagNum - 1);
  //   }
  // }
  //
  // // Convert the set to a sorted vector for sequential computation
  // std::vector<int> lagNumNeed(lagNumNeedSet.begin(), lagNumNeedSet.end());
  // std::sort(lagNumNeed.begin(), lagNumNeed.end());
  //
  // // Compute lagged neighbor results for each lagNum in the sorted sequence
  // for (int lagNum : lagNumNeed) {
  //   laggedResultsMap[lagNum] = CppLaggedNeighbor4Lattice(nb, lagNum);
  // }

  // Compute embeddings for each lag number
  for (int lagNum = startLagNum; lagNum <= endLagNum; lagNum += step) {
    // Get the lagged neighbor results for the current lagNum
    std::vector<std::vector<int>> laggedResults = laggedResultsMap[lagNum];

    // Remove duplicates with previous lagNum (if lagNum > 0)
    if (lagNum > 0) {
      std::vector<std::vector<int>> prevLaggedResults = laggedResultsMap[lagNum - 1];
      for (int i = 0; i < n; ++i) {
        // Convert previous lagged results to a set for fast lookup
        std::unordered_set<int> prevSet(prevLaggedResults[i].begin(), prevLaggedResults[i].end());

        // Remove duplicates from current lagged results
        std::vector<int> newIndices;
        for (int index : laggedResults[i]) {
          if (prevSet.find(index) == prevSet.end()) {
            newIndices.push_back(index);
          }
        }

        // If the new indices are empty, set it to a special value (e.g., std::numeric_limits<int>::min())
        if (newIndices.empty()) {
          newIndices.push_back(std::numeric_limits<int>::min());
        }

        // Update the lagged results
        laggedResults[i] = newIndices;
      }
    }

    // Compute the mean of neighbor values for each spatial unit
    for (size_t l = 0; l < laggedResults.size(); ++l) {
      std::vector<int> neighbors = laggedResults[l];

      // If the neighbors are empty or contain only the special value, leave the embedding as NaN
      if (neighbors.empty() || (neighbors.size() == 1 && neighbors[0] == std::numeric_limits<int>::min())) {
        continue;
      }

      // Compute the mean of neighbor values
      double sum = 0.0;
      int validCount = 0;

      // Loop through the neighbors to accumulate valid (non-NaN) values
      for (int idx : neighbors) {
        // Check if vec[idx] is NaN, and skip if true
        if (!std::isnan(vec[idx])) {
          sum += vec[idx];
          ++validCount;  // Increment valid count for non-NaN values
        }
      }

      // If there are valid neighbors, calculate the mean; otherwise, leave the embedding as NaN
      if (validCount > 0) {
        xEmbedings[l][(lagNum - startLagNum) / step] = sum / validCount;
      }
    }
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns;

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
      validColumns.push_back(col);
    }
  }

  // If no columns are removed, return the original xEmbedings
  if (validColumns.size() == xEmbedings[0].size()) {
    return xEmbedings;
  } else {
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

/**
 * Generates k-nearest neighbors for each spatial unit from one spatial lattice vector.
 *
 * @param vec: A vector of double values for which neighbors are to be found.
 * @param nb: A nested vector where each sub-vector contains the indices of neighbors for the corresponding element in `vec`.
 * @param k: The number of nearest neighbors to find for each element in `vec`.
 *
 * @return: A nested vector where each sub-vector contains the indices of the k-nearest neighbors for the corresponding element in `vec`.
 *
 * The function works by iterating through each element in `vec`, expanding its neighborhood by considering 1st, 2nd, ..., nth order neighbors from `nb`,
 * removing duplicates, and then selecting the k neighbors with the smallest absolute differences in value from the original element.
 */
std::vector<std::vector<int>> GenLatticeNeighbors(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    size_t k) {

  // Initialize the result vector with empty vectors
  std::vector<std::vector<int>> result(vec.size());

  // Iterate through each element in the input vector
  for (size_t i = 0; i < vec.size(); ++i) {
    // Use a set to store unique neighbor indices
    std::unordered_set<int> uniqueNeighbors;

    // Start with the direct neighbors from nb[i]
    for (int neighborIdx : nb[i]) {
      uniqueNeighbors.insert(neighborIdx);
    }

    // If the number of unique neighbors is less than k, expand the neighborhood
    if (uniqueNeighbors.size() < k) {
      // Use a queue to manage the current level of neighbors
      std::queue<int> neighborQueue;
      for (int neighborIdx : nb[i]) {
        neighborQueue.push(neighborIdx);
      }

      // Continue expanding until we have at least k unique neighbors
      while (!neighborQueue.empty() && uniqueNeighbors.size() < k) {
        int currentIdx = neighborQueue.front();
        neighborQueue.pop();

        // Add the neighbors of the current index to the queue and uniqueNeighbors set
        for (int nextNeighborIdx : nb[currentIdx]) {
          if (uniqueNeighbors.find(nextNeighborIdx) == uniqueNeighbors.end()) {
            uniqueNeighbors.insert(nextNeighborIdx);
            neighborQueue.push(nextNeighborIdx);
          }
        }
      }
    }

    // Convert the set to a vector for sorting
    std::vector<int> neighbors(uniqueNeighbors.begin(), uniqueNeighbors.end());

    // Sort the neighbors based on the absolute difference in value from the original element
    std::sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
      return std::abs(vec[a] - vec[i]) < std::abs(vec[b] - vec[i]);
    });

    // Select the top k neighbors
    if (neighbors.size() > k) {
      neighbors.resize(k);
    }

    // Store the result for the current element
    result[i] = neighbors;
  }

  return result;
}

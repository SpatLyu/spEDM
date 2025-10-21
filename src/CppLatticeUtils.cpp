#include <iostream>
#include <stdexcept>
#include <vector>
#include <queue> // for std::queue
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <unordered_map> // for std::unordered_map
#include <limits> // for std::numeric_limits
#include <cmath> // For std::isnan
#include <string>
#include "CppStats.h"

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
std::vector<std::vector<double>> CppLaggedVal4Lattice(const std::vector<double>& vec,
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
 *   vec   - A vector of values, one for each spatial unit, to be embedded.
 *   nb    - A 2D matrix where each row represents the neighbors of a spatial unit.
 *   E     - The embedding dimension, specifying how many lags to consider in the embeddings.
 *   tau   - The spatial lag step for constructing lagged state-space vectors.
 *   style - Embedding style selector:
 *             - style = 0: embedding includes current state as the first dimension.
 *             - style = 1: embedding excludes current state.
 *
 * Returns:
 *   A 2D vector representing the embeddings for each spatial unit, where each spatial unit has a row in the matrix
 *   corresponding to its embedding values for each lag number. If no valid embedding columns remain after removing
 *   columns containing only NaN values, a filtered matrix is returned.
 *
 * Note:
 *   When tau = 0, lagged variables are calculated for lag steps 0, 1, ..., E-1.
 *   When tau > 0 and style = 0, lagged variables are calculated for lag steps 0, tau, 2*tau, ..., (E-1)*tau.
 *   When tau > 0 and style != 0, lagged variables are calculated for lag steps tau, 2*tau, ..., E*tau.
 */
std::vector<std::vector<double>> GenLatticeEmbeddings(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int E,
    int tau,
    int style = 1)
{
  // Get the number of nodes
  int n = vec.size();

  // Initialize the embeddings matrix with NaN values
  std::vector<std::vector<double>> xEmbedings(n, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));

  // Precompute lagged neighbor results for all required lagNum values
  std::unordered_map<int, std::vector<std::vector<int>>> laggedResultsMap;

  // Determine the range of lagNum values and lag step based on tau and style
  int startLagNum = (style == 0) ? 0 : ( (tau == 0) ? 0 : tau );

  int endLagNum = (tau == 0)
                  ? (E - 1)
                  : (style == 0 ? (E - 1) * tau : E * tau);

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
 * Constructs composite lattice embeddings by collecting raw neighbor values
 * (instead of averaging them) for each embedding lag step.
 *
 * This function extends the traditional lattice embedding generator by
 * returning all neighbor values associated with each lag step, organized as a
 * three-dimensional vector:
 *
 *   embeddings[lag_index][unit_index][neighbor_index]
 *
 * Each inner 2D matrix corresponds to one lag (embedding dimension),
 * where rows represent spatial units and columns represent the neighbor
 * values at that lag. If the number of neighbors varies across spatial units,
 * the rows are padded with NaN values to match the maximum neighbor count
 * at that lag step.
 *
 * Parameters:
 *   vec   - A vector of values for each spatial unit.
 *   nb    - A 2D integer matrix where each row lists the neighbors of a unit.
 *   E     - Embedding dimension (number of lag steps).
 *   tau   - Spatial lag step between embedding dimensions.
 *   style - Embedding style:
 *             - 0: includes the current state as the first lag (default behavior).
 *             - 1: excludes the current state.
 *
 * Returns:
 *   A 3D vector (std::vector<std::vector<std::vector<double>>>) representing
 *   the composite lattice embeddings for all lag steps. The outermost dimension
 *   corresponds to lag steps (1 to E), the middle to spatial units, and the
 *   innermost to neighbor values (NaN-padded where needed).
 *
 * Note:
 *   - When tau = 0, lag steps are computed for 0, 1, ..., E-1.
 *   - When tau > 0 and style = 0, lag steps are computed for 0, tau, 2*tau, ..., (E-1)*tau.
 *   - When tau > 0 and style != 0, lag steps are computed for tau, 2*tau, ..., E*tau.
 */
std::vector<std::vector<std::vector<double>>> GenLatticeEmbeddingsCom(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int E,
    int tau,
    int style = 1)
{
  int n = vec.size();
  std::unordered_map<int, std::vector<std::vector<int>>> laggedResultsMap;

  // Determine lagNum range and stepping scheme
  int startLagNum = (style == 0) ? 0 : ((tau == 0) ? 0 : tau);
  int endLagNum = (tau == 0)
    ? (E - 1)
    : (style == 0 ? (E - 1) * tau : E * tau);
  int step = (tau == 0) ? 1 : tau;

  // Step 1: Build lagged neighbor index maps
  for (int lagNum = 0; lagNum <= endLagNum; ++lagNum) {
    if (lagNum == 0) {
      std::vector<std::vector<int>> result_temp;
      for (size_t i = 0; i < nb.size(); ++i)
        result_temp.push_back({static_cast<int>(i)});
      laggedResultsMap[lagNum] = result_temp;
    } else {
      std::vector<std::vector<int>> prevResult = laggedResultsMap[lagNum - 1];
      std::vector<std::vector<int>> currentResult;
      currentResult.reserve(n);

      for (int i = 0; i < n; ++i) {
        std::unordered_set<int> mergedSet;

        for (int elem : prevResult[i])
          if (elem != std::numeric_limits<int>::min())
            mergedSet.insert(elem);

          std::unordered_set<int> newElements;
          for (int j : prevResult[i]) {
            if (j == std::numeric_limits<int>::min() || j < 0 || j >= n)
              continue;
            for (int k : nb[j])
              newElements.insert(k);
          }

          mergedSet.insert(newElements.begin(), newElements.end());

          std::vector<int> vecRes(mergedSet.begin(), mergedSet.end());
          std::sort(vecRes.begin(), vecRes.end());
          vecRes.erase(std::unique(vecRes.begin(), vecRes.end()), vecRes.end());

          if (vecRes.empty())
            vecRes.push_back(std::numeric_limits<int>::min());

          currentResult.push_back(vecRes);
      }
      laggedResultsMap[lagNum] = currentResult;
    }
  }

  // Step 2: Compute embeddings for each lagNum
  std::vector<std::vector<std::vector<double>>> allLagEmbeddings;

  for (int lagNum = startLagNum; lagNum <= endLagNum; lagNum += step) {
    std::vector<std::vector<int>> laggedResults = laggedResultsMap[lagNum];

    // Deduplicate with previous lagNum
    if (lagNum > 0) {
      std::vector<std::vector<int>> prevLaggedResults = laggedResultsMap[lagNum - 1];
      for (int i = 0; i < n; ++i) {
        std::unordered_set<int> prevSet(prevLaggedResults[i].begin(),
                                        prevLaggedResults[i].end());
        std::vector<int> newIndices;
        for (int idx : laggedResults[i])
          if (prevSet.find(idx) == prevSet.end())
            newIndices.push_back(idx);

          if (newIndices.empty())
            newIndices.push_back(std::numeric_limits<int>::min());

          laggedResults[i] = newIndices;
      }
    }

    // Step 3: Determine max neighbor count across all spatial units for padding
    size_t maxNeighborCount = 0;
    for (const auto& row : laggedResults)
      maxNeighborCount = std::max(maxNeighborCount, row.size());

    // Step 4: Construct the 2D embedding matrix for this lag
    std::vector<std::vector<double>> lagMatrix(n,
                                               std::vector<double>(maxNeighborCount,
                                                                   std::numeric_limits<double>::quiet_NaN()));

    for (int i = 0; i < n; ++i) {
      const std::vector<int>& neighbors = laggedResults[i];
      for (size_t j = 0; j < neighbors.size(); ++j) {
        int idx = neighbors[j];
        if (idx >= 0 && idx < n && !std::isnan(vec[idx])) {
          lagMatrix[i][j] = vec[idx];
        } else {
          lagMatrix[i][j] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }

    allLagEmbeddings.push_back(lagMatrix);
  }

  return allLagEmbeddings;
}

/**
 * @brief Generate a list of k nearest neighbors for each spatial location based on lattice connectivity.
 *
 * This function constructs neighborhood information for each element in a spatial process
 * using both direct connectivity and value similarity. It ensures that each location has
 * at least k unique neighbors by expanding through its neighbors' neighbors recursively,
 * if necessary. All neighbors must be indices present in the provided `lib` vector.
 *
 * The procedure consists of:
 * 1. Starting with directly connected neighbors from `nb` that are also in `lib`.
 * 2. If fewer than k unique neighbors are found, iteratively expand the neighborhood using
 *    a breadth-first search (BFS) on the adjacency list (only considering nodes in `lib`).
 * 3. Among all collected neighbors, the function selects the k most similar ones in terms of
 *    absolute value difference from the center location.
 *
 * @param vec A vector of values representing the spatial process (used for sorting by similarity).
 * @param nb A list of adjacency lists where `nb[i]` gives the direct neighbors of location i.
 * @param lib A vector of indices representing valid neighbors to consider for all locations.
 * @param k The desired number of neighbors for each location.
 *
 * @return A vector of vectors, where each subvector contains the indices of the k nearest neighbors
 *         for each location, based on lattice structure and value similarity.
 */
std::vector<std::vector<int>> GenLatticeNeighbors(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& lib,
    size_t k) {

  // Preconvert lib into a set for fast lookup
  std::unordered_set<int> libSet(lib.begin(), lib.end());

  // // Check whether indices in lib are valid (recommended for robustness)
  // for (int idx : lib) {
  //   if (idx < 0 || idx >= static_cast<int>(vec.size())) {
  //     throw std::invalid_argument("Invalid index " + std::to_string(idx) + " found in 'lib'");
  //   }
  // }

  std::vector<std::vector<int>> result(vec.size());

  for (size_t i = 0; i < vec.size(); ++i) {
    std::unordered_set<int> uniqueNeighbors;

    // Initial stage: collect directly connected neighbors that exist in lib
    for (int neighborIdx : nb[i]) {
      if (libSet.count(neighborIdx)) {
        uniqueNeighbors.insert(neighborIdx);
      }
    }

    // If direct neighbors are not enough, expand using BFS (only nodes in lib)
    if (uniqueNeighbors.size() < k) {
      std::queue<int> neighborQueue;

      // Initialize the queue with valid direct neighbors
      for (int neighborIdx : nb[i]) {
        if (libSet.count(neighborIdx)) {
          neighborQueue.push(neighborIdx);
        }
      }

      // Expand neighbors using BFS until we reach k or cannot expand further
      while (!neighborQueue.empty() && uniqueNeighbors.size() < k) {
        int currentIdx = neighborQueue.front();
        neighborQueue.pop();

        // Traverse neighbors of current node and add new valid ones
        for (int nextNeighborIdx : nb[currentIdx]) {
          if (libSet.count(nextNeighborIdx) &&
              uniqueNeighbors.find(nextNeighborIdx) == uniqueNeighbors.end()) {
            uniqueNeighbors.insert(nextNeighborIdx);
            neighborQueue.push(nextNeighborIdx);
          }
        }
      }
    }

    // // Check whether enough neighbors were found
    // if (uniqueNeighbors.size() < k) {
    //   throw std::runtime_error("Location " + std::to_string(i) +
    //                            " cannot find enough (" + std::to_string(k) +
    //                            ") valid neighbors from the provided 'lib' set");
    // }

    // Convert the set to a vector and sort by value similarity
    std::vector<int> neighbors(uniqueNeighbors.begin(), uniqueNeighbors.end());
    std::sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
      return std::abs(vec[a] - vec[i]) < std::abs(vec[b] - vec[i]);
    });

    // Keep only the top-k most similar neighbors
    if (neighbors.size() > k) {
      neighbors.resize(k);
    }

    result[i] = neighbors;
  }

  return result;
}

/**
 * @brief Generate symbolization values for a spatial cross-sectional series using a lattice-based
 *        neighborhood approach, based on the method described by Herrera et al. (2016).
 *
 * This function implements a symbolic transformation of a univariate spatial process,
 * where each spatial location is associated with a value from the original series and
 * its surrounding neighborhood. The symbolization is based on comparing local median-based
 * indicators within a defined spatial neighborhood.
 *
 * The procedure follows three main steps:
 * 1. Compute the median of the input series `vec` using only the indices specified in `lib`.
 * 2. For each location in `vec`, define a binary indicator (`tau_s`) which is 1 if the value
 *    at that location is greater than or equal to the `lib`-based median, and 0 otherwise.
 * 3. For each location in `pred`, compare its indicator with those of its k nearest neighbors.
 *    The final symbolic value is the count of neighbors that share the same indicator value.
 *
 * @param vec A vector of double values representing the spatial process.
 * @param nb A nested vector containing neighborhood information (e.g., lattice connectivity).
 * @param lib A vector of indices representing valid neighbors to consider for computing the median and selecting neighbors.
 * @param pred A vector of indices specifying which elements to compute the symbolization for.
 * @param k The number of nearest neighbors to consider for each location.
 *
 * @return A vector of symbolic values (as double) for each spatial location specified in `pred`.
 */
std::vector<double> GenLatticeSymbolization(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    size_t k) {
  // Initialize the result vector with the same size as pred
  std::vector<double> result(pred.size());

  // Generate neighbors for the elements in pred
  std::vector<std::vector<int>> neighbors = GenLatticeNeighbors(vec, nb, lib, k);

  // Compute the median using only values at lib indices
  std::vector<double> lib_vals(lib.size());
  // No need to filter no-nan value
  // for (int idx : lib) {
  //   if (!std::isnan(vec[idx])) {
  //     lib_vals.push_back(vec[idx]);
  //   }
  // }
  for (size_t i = 0; i < lib.size(); ++i){
    lib_vals[i] = vec[lib[i]];
  }
  double vec_me = CppMedian(lib_vals, true);

  // Define tau_s for all positions in vec
  std::vector<double> tau_s(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    tau_s[i] = (vec[i] >= vec_me) ? 1.0 : 0.0;
  }

  // For each location in pred, compute fs
  for (size_t s = 0; s < pred.size(); ++s) {
    int currentIndex = pred[s];
    const std::vector<int>& local_neighbors = neighbors[currentIndex];
    double taus = tau_s[currentIndex];

    // Count how many neighbors share the same binary indicator
    double fs = 0.0;
    for (size_t i = 0; i < local_neighbors.size(); ++i) {
      if (tau_s[local_neighbors[i]] == taus) fs += 1.0;
    }
    result[s] = fs;
  }

  return result;
}

/**
 * @brief Divide a spatial lattice into connected blocks of approximately equal size.
 *
 * This function partitions a spatial domain represented by an adjacency list (neighbor structure)
 * into `b` spatially contiguous blocks. It ensures that each block is connected and handles isolated
 * units by merging them into the smallest neighboring block.
 *
 * @param nb A vector of vectors representing the adjacency list (i.e., neighboring indices)
 *           for each spatial unit; `nb[i]` contains the indices of neighbors of unit `i`.
 * @param b  The number of blocks to divide the lattice into.
 *
 * @return A vector of integers of length `N` where each entry corresponds to the assigned block label
 *         (ranging from 0 to b-1) of the spatial unit at that index.
 */
std::vector<int> CppDivideLattice(const std::vector<std::vector<int>>& nb, int b) {
  int N = static_cast<int>(nb.size());
  int base_size = N / b;
  int surplus = N % b;

  std::vector<bool> visited(N, false);
  std::vector<int> labels(N, -1);
  int current_block = 0;

  // Step 1: Divide the lattice into `b` blocks using BFS
  while (current_block < b) {
    int target_size = (current_block == b - 1) ? (base_size + surplus) : base_size;

    // Find the next unvisited starting point with the highest degree
    int start = -1;
    int max_degree = -1;
    for (int i = 0; i < N; ++i) {
      if (!visited[i] && static_cast<int>(nb[i].size()) > max_degree) {
        start = i;
        max_degree = nb[i].size();
      }
    }
    if (start == -1) break; // no more unvisited nodes

    std::queue<int> q;
    std::vector<int> block_members;

    q.push(start);
    visited[start] = true;

    // Perform BFS to collect connected nodes for the current block
    while (!q.empty() && static_cast<int>(block_members.size()) < target_size) {
      int current = q.front();
      q.pop();

      block_members.push_back(current);

      for (int neighbor : nb[current]) {
        if (!visited[neighbor]) {
          visited[neighbor] = true;
          q.push(neighbor);
        }
      }
    }

    // If not enough neighbors collected, prioritize connected unvisited nodes
    for (int i = 0; static_cast<int>(block_members.size()) < target_size && i < N; ++i) {
      if (!visited[i]) {
        // Check if the node is connected to any block member
        bool is_connected = false;
        for (int member : block_members) {
          if (std::find(nb[member].begin(), nb[member].end(), i) != nb[member].end()) {
            is_connected = true;
            break;
          }
        }
        if (is_connected) {
          visited[i] = true;
          block_members.push_back(i);
        }
      }
    }

    // Assign label to all members of the current block
    for (int idx : block_members) {
      labels[idx] = current_block;
    }

    ++current_block;
  }

  // Step 2: Check for isolated units and merge them into the smallest neighboring block
  for (int i = 0; i < N; ++i) {
    bool is_isolated = true;
    for (int neighbor : nb[i]) {
      if (labels[neighbor] == labels[i]) {
        is_isolated = false;
        break;
      }
    }
    if (is_isolated) {
      // Find the smallest block among neighbors
      int smallest_block = b; // Initialize with an invalid value
      for (int neighbor : nb[i]) {
        if (labels[neighbor] != -1 && labels[neighbor] < smallest_block) {
          smallest_block = labels[neighbor];
        }
      }
      if (smallest_block != b) {
        labels[i] = smallest_block; // Merge into the smallest neighboring block
      }
    }
  }

  return labels;
}

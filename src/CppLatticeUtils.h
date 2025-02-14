#ifndef CppLatticeUtils_H
#define CppLatticeUtils_H

#include <iostream>
#include <vector>
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
#include <unordered_map> // for std::unordered_map
#include <limits> // for std::numeric_limits
#include <cmath> // For std::isnan

/**
 * Computes the lagged neighbors for a lattice structure up to a given lag number.
 *
 * For lagNum=0, returns each node's index as its own neighbor.
 * For lagNum>=1, recursively expands neighbors by looking up each previous level's
 * neighbors, combines all results up to lagNum, and deduplicates. Empty results are
 * filled with std::numeric_limits<int>::min().
 *
 * @note This returns an accumulated sequence of neighbor indices.
 *
 * @param spNeighbor A 2D vector where each element contains indices of immediate neighbors.
 * @param lagNum The number of lag steps to compute.
 * @return A 2D vector of lagged neighbors for each spatial unit.
 */
std::vector<std::vector<int>> CppLaggedNeighbor4Lattice(const std::vector<std::vector<int>>& spNeighbor,
                                                        int lagNum);

/**
 * @brief Computes the lagged values for a given vector based on the neighborhood and lag number.
 *
 * This function first computes the lagged neighbors for each point in the lattice using the `CppLaggedNeighbor4Lattice` function.
 * Then, it uses these indices to extract the corresponding values from the input vector `vec`.
 *
 * @param vec The input vector of double values for which lagged values are to be computed.
 * @param nb The neighborhood matrix, where each row represents the neighbors of a unit in the spatial lattice data.
 * @param lagNum The number of lags to consider.
 * @return A vector of vectors containing the lagged values for each unit in the spatial lattice data.
 */
std::vector<std::vector<double>> CppLaggedVar4Lattice(const std::vector<double>& vec,
                                                      const std::vector<std::vector<int>>& nb,
                                                      int lagNum);

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
    int tau);

#endif // CppLatticeUtils_H

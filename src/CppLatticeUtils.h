#ifndef CppLatticeUtils_H
#define CppLatticeUtils_H

#include <iostream>
#include <vector>
#include <numeric>   // for std::accumulate
#include <algorithm> // for std::sort, std::unique, std::accumulate
#include <unordered_set> // for std::unordered_set
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
 * @param spNeighbor A 2D vector where each element contains indices of immediate neighbors.
 * @param lagNum The number of lag steps to compute.
 * @return A 2D vector of lagged neighbors for each spatial unit.
 */
std::vector<std::vector<int>> CppLaggedNeighbor4Lattice(const std::vector<std::vector<int>>& spNeighbor,
                                                        int lagNum);

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

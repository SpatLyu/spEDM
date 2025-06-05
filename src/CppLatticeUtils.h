#ifndef CppLatticeUtils_H
#define CppLatticeUtils_H

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
                                                        int lagNum);

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
    size_t k);

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
    size_t k);

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
std::vector<int> CppDivideLattice(
    const std::vector<std::vector<int>>& nb,
    int b);

#endif // CppLatticeUtils_H

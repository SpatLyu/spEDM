#ifndef SpatialBlockBootstrap_H
#define SpatialBlockBootstrap_H

#include <vector>
#include <unordered_map>
#include <random>

/**
 * @brief Generate a spatial block bootstrap resample of block indices based on predefined blocks.
 *
 * This function follows the Spatial Block Bootstrap (SBB) procedure described by Carlstein (1986)
 * and Herrera et al. (2013), as used in spatial Granger causality frameworks. It samples blocks
 * of indices with replacement, preserving spatial contiguity defined by a block ID vector.
 *
 *
 * @param block     Vector of block IDs assigning each element to a contiguous block.
 *                  Elements with the same integer value belong to the same block.
 * @param seed      Random seed for reproducibility (optional).
 * @return std::vector<int> The bootstrap-resampled vector of indices based on block IDs.
 */
std::vector<int> SpatialBlockBootstrap(
    const std::vector<int>& block,
    unsigned int seed = 42
);

#endif // SpatialBlockBootstrap_H

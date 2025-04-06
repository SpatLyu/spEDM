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
) {
  size_t N = block.size();

  // Step 1: Group indices by block ID
  std::unordered_map<int, std::vector<size_t>> block_to_indices;
  for (size_t i = 0; i < N; ++i) {
    block_to_indices[block[i]].push_back(i);
  }

  // Step 2: Extract unique block IDs
  std::vector<int> block_ids;
  for (const auto& pair : block_to_indices) {
    block_ids.push_back(pair.first);
  }

  // Step 3: Random sampling of block IDs with replacement
  std::mt19937 rng(seed); // Random number generator with fixed seed
  std::uniform_int_distribution<> dist(0, block_ids.size() - 1);

  // Step 4: Generate bootstrap sample by sampling blocks
  std::vector<int> bootstrapped_indices;
  while (bootstrapped_indices.size() < N) {
    int sampled_block_id = block_ids[dist(rng)];
    const auto& indices = block_to_indices[sampled_block_id];

    // Append block indices to the bootstrap sample
    for (size_t idx : indices) {
      bootstrapped_indices.push_back(idx);
    }
  }

  // Trim to original size (in case last block pushed size beyond N)
  if (bootstrapped_indices.size() > N) {
    bootstrapped_indices.resize(N);
  }

  return bootstrapped_indices;
}

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

/**
 * Generate a spatial block bootstrap sample of indices using an external random number generator.
 *
 * This function performs block-based resampling from spatial or grouped data. The input vector `block`
 * specifies a block ID for each observation (e.g., spatial unit, group, or time block). The function:
 *
 *   1. Groups indices by block ID;
 *   2. Randomly samples block IDs with replacement using the provided `std::mt19937` RNG;
 *   3. Concatenates the indices of the selected blocks to form the bootstrap sample;
 *   4. Trims the result to ensure the same length as the original data.
 *
 * The sampling is done at the block level rather than at the individual level, preserving local structure.
 * This is particularly useful for spatial or temporal data where neighboring observations may be dependent.
 *
 * @param block  A vector of block identifiers (one per observation), e.g., spatial or temporal blocks.
 * @param rng    A reference to an externally managed random number generator (e.g., from a parallel RNG pool).
 * @return       A vector of resampled indices with the same length as the input data.
 */
std::vector<int> SpatialBlockBootstrapRNG(
    const std::vector<int>& block,
    std::mt19937_64& rng
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

  // Step 3: Sampling of block IDs with replacement using external RNG
  std::uniform_int_distribution<> dist(0, block_ids.size() - 1);

  // Step 4: Generate bootstrap sample
  std::vector<int> bootstrapped_indices;
  while (bootstrapped_indices.size() < N) {
    int sampled_block_id = block_ids[dist(rng)];
    const auto& indices = block_to_indices[sampled_block_id];
    bootstrapped_indices.insert(bootstrapped_indices.end(), indices.begin(), indices.end());
  }

  // Trim to original size if needed
  if (bootstrapped_indices.size() > N) {
    bootstrapped_indices.resize(N);
  }

  return bootstrapped_indices;
}

#include <vector>
#include <cmath>
#include <limits>
#include "CppStats.h"

/*
 * Compute the False Nearest Neighbors (FNN) ratio for spatial cross-sectional data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - Rtol: Relative threshold (default 10.0). If the change in the added dimension is
 *         large relative to the E1 distance, the neighbor is considered false.
 * - Atol: Absolute threshold (default 2.0). If the added dimension changes too much
 *         absolutely, the neighbor is also considered false.
 * - L1norm: Whether to use Manhattan (L1) distance instead of Euclidean (L2).
 *
 * Returns:
 * - A double value indicating the proportion of false nearest neighbors (0–1).
 *   If no valid pairs are found, returns NaN.
 */
double CppSingleFNN(const std::vector<std::vector<double>>& embedding,
                    size_t E1,
                    size_t E2,
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false) {
  size_t N = embedding.size();

  int false_count = 0;
  int total = 0;

  if (embedding.empty() || embedding[0].size() < E2) {
    return std::numeric_limits<double>::quiet_NaN();  // Invalid dimensions
  }

  for (size_t i = 0; i < N; ++i) {
    if (checkOneDimVectorNotNanNum(embedding[i]) == 0) {
      continue;  // Skip rows with all NaNs
    }

    // Extract E1-dimensional embedding for unit i
    std::vector<double> xi_E1(embedding[i].begin(), embedding[i].begin() + E1);

    double min_dist = std::numeric_limits<double>::max();
    size_t nn_idx = N;  // invalid index placeholder

    // Find nearest neighbor of i in E1-dimensional space
    for (size_t j = 0; j < N; ++j) {
      if (i == j || checkOneDimVectorNotNanNum(embedding[j]) == 0) continue;

      std::vector<double> xj_E1(embedding[j].begin(), embedding[j].begin() + E1);
      double dist = CppDistance(xi_E1, xj_E1, L1norm, true);  // true: skip NaNs

      if (dist < min_dist) {
        min_dist = dist;
        nn_idx = j;
      }
    }

    if (nn_idx == N || min_dist == 0.0) continue;  // skip degenerate cases

    // Compare E2-th coordinate difference (new dimension)
    double diff = std::abs(embedding[i][E2 - 1] - embedding[nn_idx][E2 - 1]);
    double ratio = diff / min_dist;

    if (ratio > Rtol || diff > Atol) {
      ++false_count;
    }
    ++total;
  }

  return total > 0 ? static_cast<double>(false_count) / total
  : std::numeric_limits<double>::quiet_NaN();
}

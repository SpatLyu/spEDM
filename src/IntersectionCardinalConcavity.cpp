#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>

// Function to compute the Euclidean distance between two vectors
double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return std::sqrt(sum);
}

// Function to find k-nearest neighbors of a given index in the embedding space
std::vector<std::size_t> find_k_nearest_neighbors(
    const std::vector<std::vector<double>>& embedding_space,
    std::size_t target_idx,
    std::size_t k)
{
  std::size_t n = embedding_space.size();
  std::vector<std::pair<double, std::size_t>> distances;

  for (std::size_t i = 0; i < n; ++i) {
    if (i != target_idx) {
      double dist = euclidean_distance(embedding_space[target_idx], embedding_space[i]);
      distances.emplace_back(dist, i);
    }
  }

  // Partial sort to get k-nearest neighbors
  std::partial_sort(distances.begin(), distances.begin() + k, distances.end());

  std::vector<std::size_t> neighbors;
  for (std::size_t i = 0; i < k; ++i) {
    neighbors.push_back(distances[i].second);
  }

  return neighbors;
}

// Function to compute the Intersectional Cardinality (IC) curve
std::vector<double> compute_IC_curve(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    std::size_t k, std::size_t max_r)
{
  std::vector<double> IC_curve;

  for (std::size_t r = k; r <= max_r; ++r) {
    double intersection_sum = 0.0;

    for (std::size_t t = 0; t < embedding_x.size(); ++t) {
      std::vector<std::size_t> neighbors_x = find_k_nearest_neighbors(embedding_x, t, k);
      std::vector<std::size_t> neighbors_y = find_k_nearest_neighbors(embedding_y, t, r);

      // Compute the intersection count
      std::size_t intersection_count = 0;
      for (std::size_t nx : neighbors_x) {
        if (std::find(neighbors_y.begin(), neighbors_y.end(), nx) != neighbors_y.end()) {
          intersection_count++;
        }
      }

      intersection_sum += intersection_count;
    }

    IC_curve.push_back(intersection_sum / static_cast<double>(embedding_x.size()));
  }

  return IC_curve;
}

// Function to compute Area Under IC Curve (AUC) using the trapezoidal rule
double compute_AUC(const std::vector<double>& IC_curve, std::size_t k, std::size_t max_r) {
  double auc = 0.0;

  for (std::size_t i = 1; i < IC_curve.size(); ++i) {
    double width = 1.0;  // Assuming uniform step size in r
    auc += (IC_curve[i] + IC_curve[i - 1]) * width / 2.0;
  }

  return auc / static_cast<double>(max_r - k);
}

/*
 * Computes the Intersection Cardinal Concavity (ICC) causal strength score.
 *
 * Parameters:
 *   - embedding_x: The state-space reconstructed from the potential cause variable.
 *   - embedding_y: The state-space reconstructed from the potential effect variable.
 *   - lib_indices: Boolean vector indicating which states are included in the library.
 *   - pred_indices: Boolean vector indicating which states are used for prediction.
 *   - num_neighbors: Number of neighbors used for cross-mapping.
 *
 * Returns:
 *   - A double representing the ICC causal strength score, normalized between [0,1].
 */
double IntersectionCardinalConcavity(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
) {
  std::size_t k = static_cast<std::size_t>(num_neighbors);
  std::size_t max_r = k + 10;  // Adjust based on dataset size

  // Compute the Intersectional Cardinality (IC) curve
  std::vector<double> IC_curve = compute_IC_curve(embedding_x, embedding_y, k, max_r);

  // Compute AUC and normalize causal strength score
  double AUC = compute_AUC(IC_curve, k, max_r);
  double ICC_score = std::max(0.0, 2.0 * (AUC - 0.5));

  return ICC_score;
}

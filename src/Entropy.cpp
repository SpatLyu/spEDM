#include <cmath>
#include <vector>
#include "CppStats.h"

double CppEntropy(const std::vector<double>& vec,
                  size_t k, double base = 10,
                  bool L1norm = false, bool NA_rm = false) {
  std::vector<double> distances = CppKNearestDistance(vec, k, L1norm, NA_rm);
  size_t n = vec.size();

  double sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    sum += CppLog(2 * distances[i], base);  // Apply logarithm transformation
  }
  sum /= n;

  // Compute entropy using CppDigamma function
  double E = CppDigamma(n) - CppDigamma(k) + sum + CppLog(1.0, base);
  return E;
}

double CppJoinEntropy(const std::vector<std::vector<double>>& mat,
                      size_t k, double base = 10,
                      bool L1norm = false, bool NA_rm = false) {
  size_t nrow = mat.size();
  size_t ncol = mat[0].size();

  std::vector<double> distances(nrow);
  std::vector<std::vector<double>> mat_dist = CppMatDistance(mat, L1norm, NA_rm);

  for (size_t i = 0; i < nrow; ++i) {
    // Create a vector to store the distances for the current row, filtering out NaN values if NA_rm is true
    std::vector<double> dist_n;

    if (NA_rm) {
      for (double val : mat_dist[i]) {
        if (!std::isnan(val)) {
          dist_n.push_back(val);  // Only include non-NaN values
        }
      }
    } else {
      dist_n = mat_dist[i];  // Include all values if NA_rm is false
    }

    // Use nth_element to partially sort the distances up to the k-th element
    // This is more efficient than fully sorting the entire vector.
    if (k < dist_n.size()) {
      std::nth_element(dist_n.begin(), dist_n.begin() + k, dist_n.end());
      distances[i] = dist_n[k];  // `k` is 0-indexed, so this is the (k+1)-th smallest distance
    } else {
      distances[i] = *std::max_element(dist_n.begin(), dist_n.end());  // Handle case where k is out of bounds
    }
  }

  double sum = 0.0;
  for (size_t i = 0; i < nrow; i++) {
    sum += CppLog(2 * distances[i], base);  // Apply logarithm transformation
  }
  sum = sum * ncol / nrow;

  // Compute joint entropy using CppDigamma function
  double E = CppDigamma(nrow) - CppDigamma(k) + sum;
  return E;
}

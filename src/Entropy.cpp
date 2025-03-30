#include <cmath>
#include <vector>
#include "CppStats.h"

// Alias for convenience
using VecD = std::vector<double>;
using MatD = std::vector<std::vector<double>>;


double CppEntropy(const std::vector<double>& vec,
                  int k, double base = 10,
                  bool L1norm = false, bool NA_rm = false) {
  std::vector<double> distances = CppKNearestDistance(vec, k, L1norm, NA_rm);
  size_t n = vec.size();

  double sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    sum += CppLog(2 * distances[i], base);  // Apply logarithm transformation
  }
  sum /= n;

  // Compute entropy using digamma function
  double E = std::tgamma(n) - std::tgamma(k) + sum + CppLog(1.0, base);
  return E;
}

/**
 * @brief Computes the joint entropy of a dataset represented as a matrix.
 *
 * @param M A matrix where each row represents a data point.
 * @param k The number of nearest neighbors to consider.
 * @param log_type The logarithm type (e.g., "ln", "log2", "log10").
 * @return The computed joint entropy value.
 */
double CppJoinEntropy(const MatD& M, int k, const std::string& log_type) {
  double E = 0.0;
  unsigned N = M.size();
  double sum = 0.0;

  unsigned d = M[0].size();  // Dimensionality of the data

  VecD distances = kNearest(M, k);  // Compute k-nearest neighbor distances

  for (unsigned i = 0; i < N; i++) {
    sum += myLOG(2 * distances[i], log_type);  // Apply logarithm transformation
  }

  sum = sum * d / N;

  // Compute joint entropy using digamma function
  E = std::tgamma(N) - std::tgamma(k) + sum;
  return E;
}

#include <cmath>
#include <vector>
#include "CppStats.h"

double CppEntropy(const std::vector<double>& vec, size_t k,
                  double base = 10, bool NA_rm = false) {
  std::vector<double> distances = CppKNearestDistance(vec, k, true, NA_rm);
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

double CppJoinEntropy(const std::vector<std::vector<double>>& mat, size_t k,
                      double base = 10, bool NA_rm = false) {
  size_t nrow = mat.size();
  size_t ncol = mat[0].size();

  std::vector<double> distances(nrow);
  std::vector<std::vector<double>> mat_dist = CppMatDistance(mat, true, NA_rm);

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
      distances[i] = dist_n[k];  // (k+1)-th smallest distance (exclude itself)
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

double CppMutualInformation(const std::vector<std::vector<double>>& mat, size_t k,
                            int alg = 1, bool normalize = false, bool NA_rm = false){
  size_t nrow = mat.size();
  // size_t ncol = mat[0].size();

  std::vector<double> X(nrow);
  std::vector<double> Y(nrow);
  for (size_t i = 0; i < nrow; ++i) {
    X[i] = mat[i][0];
    Y[i] = mat[i][1];
  }

  std::vector<double> distances(nrow);
  std::vector<std::vector<double>> mat_dist = CppMatDistance(mat, true, NA_rm);

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
      distances[i] = dist_n[k];  // (k+1)-th smallest distance (exclude itself)
    } else {
      distances[i] = *std::max_element(dist_n.begin(), dist_n.end());  // Handle case where k is out of bounds
    }
  }

  double sum = 0;
  double mi = 0;
  if (alg == 1){
    std::vector<int> NX = CppNeighborsNum(X, distances, false, true, NA_rm);
    std::vector<int> NY = CppNeighborsNum(Y, distances, false, true, NA_rm);
    for (size_t i = 0; i < nrow; i ++){
      sum += CppDigamma(NX[i] + 1) + CppDigamma(NY[i] + 1);
    }
    sum /= nrow;
    mi = CppDigamma(k) + CppDigamma(nrow) - sum;
  } else {
    std::vector<double> distances_x = CppKNearestDistance(X, k, true, NA_rm);
    std::vector<double> distances_y = CppKNearestDistance(Y, k, true, NA_rm);
    std::vector<int> NX = CppNeighborsNum(X, distances_x, true, true, NA_rm);
    std::vector<int> NY = CppNeighborsNum(Y, distances_y, true, true, NA_rm);
    for (size_t i = 0; i < nrow; i++){
      sum += CppDigamma(NX[i]) + CppDigamma(NY[i]);
    }
    sum /= nrow;
    mi = CppDigamma(k) - (1.0 / k) + CppDigamma(nrow) - sum;
  }

  // Normalizing mutual information by divide it by the joint entropy
  if (normalize) {
    double jointEn = 0;
    for (double d: distances){
      jointEn += d;
    }
    jointEn *= (2.0 / distances.size());
    jointEn +=  CppDigamma(nrow) - CppDigamma(k);
    mi = mi / jointEn;
  }
  return mi;
}

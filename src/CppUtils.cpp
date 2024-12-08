#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate

// Function to calculate the mean of a vector
double CppMean(const std::vector<double>& vec) {
  return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

// Function to calculate the variance of a vector
double CppVariance(const std::vector<double>& vec) {
  double mean_val = CppMean(vec);
  double var = 0.0;
  for (const auto& value : vec) {
    var += (value - mean_val) * (value - mean_val);
  }
  return var / vec.size();
}

// Function to calculate the covariance of two vectors
double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  double mean1 = CppMean(vec1);
  double mean2 = CppMean(vec2);
  double cov = 0.0;
  for (size_t i = 0; i < vec1.size(); ++i) {
    cov += (vec1[i] - mean1) * (vec2[i] - mean2);
  }
  return cov / vec1.size();
}

// Function to calculate the Pearson correlation coefficient
double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat) {
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  double cov_yy_hat = CppCovariance(y, y_hat);
  double var_y = CppVariance(y);
  double var_y_hat = CppVariance(y_hat);

  return cov_yy_hat / std::sqrt(var_y * var_y_hat);
}

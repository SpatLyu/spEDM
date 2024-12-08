#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits

// Function to check if a value is NA
bool isNA(double value) {
  return std::isnan(value);
}

// Function to calculate the mean of a vector, ignoring NA values
double CppMean(const std::vector<double>& vec, bool NA_rm = false) {
  double sum = 0.0;
  size_t count = 0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
      ++count;
    }
  }
  return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the variance of a vector, ignoring NA values
double CppVariance(const std::vector<double>& vec, bool NA_rm = false) {
  double mean_val = CppMean(vec, NA_rm);
  double var = 0.0;
  size_t count = 0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      var += (value - mean_val) * (value - mean_val);
      ++count;
    }
  }
  return count > 1 ? var / (count - 1) : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the covariance of two vectors, ignoring NA values
double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2,
                     bool NA_rm = false) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  double mean1 = CppMean(vec1, NA_rm);
  double mean2 = CppMean(vec2, NA_rm);
  double cov = 0.0;
  size_t count = 0;
  for (size_t i = 0; i < vec1.size(); ++i) {
    if ((!NA_rm || !isNA(vec1[i])) && (!NA_rm || !isNA(vec2[i]))) {
      cov += (vec1[i] - mean1) * (vec2[i] - mean2);
      ++count;
    }
  }
  return count > 1 ? cov / (count - 1) : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the Pearson correlation coefficient, ignoring NA values
double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false) {
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  double cov_yy_hat = CppCovariance(y, y_hat, NA_rm);
  double var_y = CppVariance(y, NA_rm);
  double var_y_hat = CppVariance(y_hat, NA_rm);

  return cov_yy_hat / std::sqrt(var_y * var_y_hat);
}

// // Function to calculate the mean of a vector
// double CppMean(const std::vector<double>& vec) {
//   return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
// }
//
// // Function to calculate the variance of a vector
// double CppVariance(const std::vector<double>& vec) {
//   double mean_val = CppMean(vec);
//   double var = 0.0;
//   for (const auto& value : vec) {
//     var += (value - mean_val) * (value - mean_val);
//   }
//   return var / vec.size();
// }
//
// // Function to calculate the covariance of two vectors
// double CppCovariance(const std::vector<double>& vec1,
//                      const std::vector<double>& vec2) {
//   if (vec1.size() != vec2.size()) {
//     throw std::invalid_argument("Vectors must have the same size");
//   }
//
//   double mean1 = CppMean(vec1);
//   double mean2 = CppMean(vec2);
//   double cov = 0.0;
//   for (size_t i = 0; i < vec1.size(); ++i) {
//     cov += (vec1[i] - mean1) * (vec2[i] - mean2);
//   }
//   return cov / vec1.size();
// }
//
// // Function to calculate the Pearson correlation coefficient
// double PearsonCor(const std::vector<double>& y,
//                   const std::vector<double>& y_hat) {
//   if (y.size() != y_hat.size()) {
//     throw std::invalid_argument("Vectors must have the same size");
//   }
//
//   double cov_yy_hat = CppCovariance(y, y_hat);
//   double var_y = CppVariance(y);
//   double var_y_hat = CppVariance(y_hat);
//
//   return cov_yy_hat / std::sqrt(var_y * var_y_hat);
// }

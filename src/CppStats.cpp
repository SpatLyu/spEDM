#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits
#include "HelperFuns.h"
#include <Rcpp.h>

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

// Function to calculate the significance of a correlation coefficient
double CppSignificance(double r, int n) {
  double t = r * std::sqrt((n - 2) / (1 - r * r));
  return (1 - R::pt(t, n - 2, true, false)) * 2;
}

// Function to calculate the confidence interval for a correlation coefficient
std::vector<double> CppConfidence(double r, int n,
                                  double level = 0.05) {
  // Calculate the Fisher's z-transformation
  double z = 0.5 * std::log((1 + r) / (1 - r));

  // Calculate the standard error of z
  double ztheta = 1 / std::sqrt(n - 3);

  // Calculate the z-value for the given confidence level
  double qZ = R::qnorm(1 - level / 2, 0.0, 1.0, true, false);

  // Calculate the upper and lower bounds of the confidence interval
  double upper = z + qZ * ztheta;
  double lower = z - qZ * ztheta;

  // Convert the bounds back to correlation coefficients
  double r_upper = (std::exp(2 * upper) - 1) / (std::exp(2 * upper) + 1);
  double r_lower = (std::exp(2 * lower) - 1) / (std::exp(2 * lower) + 1);

  // Return the result as a std::vector<double>
  return {r_upper, r_lower};
}

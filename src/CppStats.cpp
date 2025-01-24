#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to check if a value is NA
bool isNA(double value) {
  return std::isnan(value);
}

// Function to check if a indice of int type is NA
bool checkIntNA(int value) {
  return value == std::numeric_limits<int>::min();
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

// Function to calculate the sum of a vector, ignoring NA values if NA_rm is true
double CppSum(const std::vector<double>& vec,
              bool NA_rm = false) {
  double sum = 0.0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
    }
  }
  return sum;
}

// Function to calculate the absolute difference between two vectors
std::vector<double> CppAbs(const std::vector<double>& vec1,
                           const std::vector<double>& vec2) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  std::vector<double> result(vec1.size());
  for (size_t i = 0; i < vec1.size(); ++i) {
    result[i] = std::abs(vec1[i] - vec2[i]);
  }
  return result;
}

// Function to normalize a vector by dividing each element by the sum of all elements
std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false) {
  double sum = CppSum(vec, NA_rm);
  if (sum == 0.0) {
    throw std::invalid_argument("Sum of vector elements is zero, cannot normalize.");
  }

  std::vector<double> normalizedVec(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!isNA(vec[i])) {
      normalizedVec[i] = vec[i] / sum;
    } else {
      normalizedVec[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  return normalizedVec;
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

// Function to compute SVD similar to R's svd()
// Input:
//   - X: A matrix represented as std::vector<std::vector<double>>
// Output:
//   - A std::vector containing three components:
//     1. d: A vector of singular values (std::vector<double>)
//     2. u: A matrix of left singular vectors (std::vector<std::vector<double>>)
//     3. v: A matrix of right singular vectors (std::vector<std::vector<double>>)
std::vector<std::vector<std::vector<double>>> CppSVD(const std::vector<std::vector<double>>& X) {
  // Convert input matrix to Armadillo matrix
  size_t m = X.size();
  size_t n = X[0].size();
  arma::mat A(m, n);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      A(i, j) = X[i][j];
    }
  }

  // Perform SVD using Armadillo
  arma::mat U; // Left singular vectors
  arma::vec S; // Singular values
  arma::mat V; // Right singular vectors
  arma::svd(U, S, V, A);

  // Convert Armadillo objects back to std::vector
  std::vector<std::vector<double>> u(m, std::vector<double>(m));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      u[i][j] = U(i, j);
    }
  }

  std::vector<double> d(S.n_elem);
  for (size_t i = 0; i < S.n_elem; ++i) {
    d[i] = S(i);
  }

  std::vector<std::vector<double>> v(V.n_rows, std::vector<double>(V.n_cols));
  for (size_t i = 0; i < V.n_rows; ++i) {
    for (size_t j = 0; j < V.n_cols; ++j) {
      v[i][j] = V(i, j);
    }
  }

  // Return as a std::vector to match R's svd() output
  return {u, {d}, v};
}

// Function to perform Linear Trend Removal
std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false) {
  if (vec.size() != xcoord.size() || vec.size() != ycoord.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Perform linear regression
  double x1_mean = CppMean(xcoord, NA_rm);
  double x2_mean = CppMean(ycoord, NA_rm);
  double y_mean = CppMean(vec, NA_rm);

  double x1_var = CppVariance(xcoord, NA_rm);
  double x2_var = CppVariance(ycoord, NA_rm);
  double x1_x2_cov = CppCovariance(xcoord, ycoord, NA_rm);

  double x1_y_cov = CppCovariance(xcoord, vec, NA_rm);
  double x2_y_cov = CppCovariance(ycoord, vec, NA_rm);

  double denom = x1_var * x2_var - x1_x2_cov * x1_x2_cov;
  if (denom == 0.0) {
    throw std::invalid_argument("Linear regression cannot be performed due to collinearity.");
  }

  double b1 = (x2_var * x1_y_cov - x1_x2_cov * x2_y_cov) / denom;
  double b2 = (x1_var * x2_y_cov - x1_x2_cov * x1_y_cov) / denom;
  double b0 = y_mean - b1 * x1_mean - b2 * x2_mean;

  // Predict vec_hat using the linear regression model
  std::vector<double> vec_hat(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    vec_hat[i] = b0 + b1 * xcoord[i] + b2 * ycoord[i];
  }

  // Calculate vec - vec_hat
  std::vector<double> result(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    result[i] = vec[i] - vec_hat[i];
  }

  return result;
}

// Function to remove linear trend using Armadillo for internal calculations
std::vector<double> ArmaLinearTrendRM(const std::vector<double>& vec,
                                      const std::vector<double>& xcoord,
                                      const std::vector<double>& ycoord,
                                      bool NA_rm = false) {
  // Check input sizes
  if (vec.size() != xcoord.size() || vec.size() != ycoord.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Create a vector to store the result
  std::vector<double> result(vec.size(), std::numeric_limits<double>::quiet_NaN()); // Initialize with NA

  // Find indices where all three vectors are not NA
  std::vector<size_t> valid_indices;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!isNA(vec[i]) && !isNA(xcoord[i]) && !isNA(ycoord[i])) {
      valid_indices.push_back(i);
    } else if (!NA_rm) {
      throw std::invalid_argument("Input contains NA values and NA_rm is false.");
    }
  }

  // If no valid data, return the result filled with NA
  if (valid_indices.empty()) {
    return result;
  }

  // Extract non-NA values
  std::vector<double> clean_vec, clean_xcoord, clean_ycoord;
  for (size_t i : valid_indices) {
    clean_vec.push_back(vec[i]);
    clean_xcoord.push_back(xcoord[i]);
    clean_ycoord.push_back(ycoord[i]);
  }

  // Convert cleaned vectors to Armadillo vectors
  arma::vec arma_vec(clean_vec);
  arma::vec arma_xcoord(clean_xcoord);
  arma::vec arma_ycoord(clean_ycoord);

  // Create design matrix for linear regression
  arma::mat X(clean_vec.size(), 3);
  X.col(0) = arma::ones<arma::vec>(clean_vec.size()); // Intercept
  X.col(1) = arma_xcoord; // x1
  X.col(2) = arma_ycoord; // x2

  // Perform linear regression using Armadillo
  arma::vec coefficients;
  if (!arma::solve(coefficients, X, arma_vec)) {
    throw std::invalid_argument("Linear regression failed due to singular matrix.");
  }

  // Predict vec_hat using the linear regression model
  arma::vec arma_vec_hat = X * coefficients;

  // Fill the result vector
  for (size_t i = 0; i < valid_indices.size(); ++i) {
    result[valid_indices[i]] = clean_vec[i] - arma_vec_hat(i);
  }

  return result;
}

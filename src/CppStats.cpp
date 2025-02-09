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

// Function to compute Mean Absolute Error (MAE) between two vectors
double CppMAE(const std::vector<double>& vec1,
              const std::vector<double>& vec2,
              bool NA_rm = false) {
  // Check if input vectors have the same size
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Initialize variables for MAE calculation
  double sum_abs_diff = 0.0; // Sum of absolute differences
  size_t valid_count = 0;    // Count of valid (non-NaN) pairs

  // Iterate through the vectors
  for (size_t i = 0; i < vec1.size(); ++i) {
    // Check if either vec1[i] or vec2[i] is NaN
    if (isNA(vec1[i]) || isNA(vec2[i])) {
      if (!NA_rm) {
        // If NA_rm is false and NaN is encountered, return NaN
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      // If both values are valid, compute absolute difference and add to sum
      sum_abs_diff += std::fabs(vec1[i] - vec2[i]);
      valid_count++;
    }
  }

  // If no valid pairs are found, return NaN
  if (valid_count == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Compute and return MAE
  return sum_abs_diff / static_cast<double>(valid_count);
}

// Function to compute Root Mean Squared Error (RMSE) between two vectors
double CppRMSE(const std::vector<double>& vec1,
               const std::vector<double>& vec2,
               bool NA_rm = false) {
  // Check if input vectors have the same size
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Initialize variables for RMSE calculation
  double sum_squared_diff = 0.0; // Sum of squared differences
  size_t valid_count = 0;        // Count of valid (non-NaN) pairs

  // Iterate through the vectors
  for (size_t i = 0; i < vec1.size(); ++i) {
    // Check if either vec1[i] or vec2[i] is NaN
    if (isNA(vec1[i]) || isNA(vec2[i])) {
      if (!NA_rm) {
        // If NA_rm is false and NaN is encountered, return NaN
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      // If both values are valid, compute squared difference and add to sum
      double diff = vec1[i] - vec2[i];
      sum_squared_diff += std::pow(diff, 2);
      valid_count++;
    }
  }

  // If no valid pairs are found, return NaN
  if (valid_count == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Compute and return RMSE
  return std::sqrt(sum_squared_diff / static_cast<double>(valid_count));
}

// Function to calculate the absolute difference between two vectors
std::vector<double> CppAbsDiff(const std::vector<double>& vec1,
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

// Function to compute distance between two vectors:
double CppDistance(const std::vector<double>& vec1,
                   const std::vector<double>& vec2,
                   bool L1norm = false,
                   bool NA_rm = false){
  // Handle NA values
  std::vector<double> clean_v1, clean_v2;
  for (size_t i = 0; i < vec1.size(); ++i) {
    bool is_na = isNA(vec1[i]) || isNA(vec2[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_v1.push_back(vec1[i]);
      clean_v2.push_back(vec2[i]);
    }
  }

  // If no valid data, return NaN
  if (clean_v1.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double dist_res = 0.0;
  if (L1norm) {
    for (std::size_t i = 0; i < clean_v1.size(); ++i) {
      dist_res += std::abs(clean_v1[i] - clean_v2[i]);
    }
  } else {
    for (std::size_t i = 0; i < clean_v1.size(); ++i) {
      dist_res += (clean_v1[i] - clean_v2[i]) * (clean_v1[i] - clean_v2[i]);
    }
    dist_res = std::sqrt(dist_res);
  }

  return dist_res;
}

// Function to compute Pearson correlation using Armadillo
double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false) {
  // Check input sizes
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Handle NA values
  std::vector<double> clean_y, clean_y_hat;
  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
    }
  }

  // If no valid data, return NaN
  if (clean_y.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Convert cleaned vectors to Armadillo vectors
  arma::vec arma_y(clean_y);
  arma::vec arma_y_hat(clean_y_hat);

  // Compute Pearson correlation using Armadillo
  double corr = arma::as_scalar(arma::cor(arma_y, arma_y_hat));

  // Ensure correlation is within valid range [-1, 1]
  if (corr < -1.0) corr = -1.0;
  if (corr > 1.0) corr = 1.0;

  return corr;
}

// // Function to calculate the Pearson correlation coefficient, ignoring NA values
// double PearsonCor(const std::vector<double>& y,
//                   const std::vector<double>& y_hat,
//                   bool NA_rm = false) {
//   // // Check input sizes
//   // if (y.size() != y_hat.size()) {
//   //   throw std::invalid_argument("Input vectors must have the same size.");
//   // }
//
//   // Handle NA values
//   std::vector<double> clean_y, clean_y_hat;
//   for (size_t i = 0; i < y.size(); ++i) {
//     bool is_na = isNA(y[i]) || isNA(y_hat[i]);
//     if (is_na) {
//       if (!NA_rm) {
//         return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
//       }
//     } else {
//       clean_y.push_back(y[i]);
//       clean_y_hat.push_back(y_hat[i]);
//     }
//   }
//
//   // If no valid data, return NaN
//   if (clean_y.empty()) {
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   double cov_yy_hat = CppCovariance(clean_y, clean_y_hat, true);
//   double var_y = CppVariance(clean_y, true);
//   double var_y_hat = CppVariance(clean_y_hat, true);
//
//   // If any of the values is NaN, return NaN
//   if (isNA(cov_yy_hat) || isNA(var_y) || isNA(var_y_hat)) {
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   // Check if variances are zero
//   if (var_y == 0.0 || var_y_hat == 0.0) {
//     return std::numeric_limits<double>::quiet_NaN(); // Return NaN if variance is zero
//   }
//
//   // Calculate Pearson correlation coefficient
//   double corr = cov_yy_hat / std::sqrt(var_y * var_y_hat);
//
//   // Ensure correlation is within valid range [-1, 1]
//   if (corr < -1.0) corr = -1.0;
//   if (corr > 1.0) corr = 1.0;
//
//   return corr;
// }

// Function to compute Partial Correlation using Armadillo
// y: Dependent variable vector
// y_hat: Predicted variable vector
// controls: Matrix of control variables (**each row represents a control variable**)
// NA_rm: Boolean flag to indicate whether to remove NA values
// linear: Boolean flag to indicate whether to calculate the partial correlation coefficient using linear regression or correlation matrix
// Returns: Partial correlation between y and y_hat after controlling for the variables in controls
double PartialCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  const std::vector<std::vector<double>>& controls,
                  bool NA_rm = false,
                  bool linear = false) {
  // Check input sizes
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors y and y_hat must have the same size.");
  }
  if (!controls.empty() && controls[0].size() != y.size()) {
    throw std::invalid_argument("Control variables must have the same number of observations as y and y_hat.");
  }

  // Handle NA values
  std::vector<double> clean_y, clean_y_hat;
  std::vector<std::vector<double>> clean_controls(controls.size());

  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    for (const auto& control : controls) {
      if (isNA(control[i])) {
        is_na = true;
        break;
      }
    }
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
      for (size_t j = 0; j < controls.size(); ++j) {
        clean_controls[j].push_back(controls[j][i]);
      }
    }
  }

  // If no valid data, return NaN
  if (clean_y.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double partial_corr;
  if (linear){
    // Convert cleaned vectors to Armadillo vectors/matrices
    arma::vec arma_y(clean_y);
    arma::vec arma_y_hat(clean_y_hat);
    arma::mat arma_controls(clean_y.size(), controls.size());
    for (size_t i = 0; i < controls.size(); ++i) {
      arma_controls.col(i) = arma::vec(clean_controls[i]);
    }

    // Compute residuals of y and y_hat after regressing on controls
    arma::vec residuals_y = arma_y - arma_controls * arma::solve(arma_controls, arma_y);
    arma::vec residuals_y_hat = arma_y_hat - arma_controls * arma::solve(arma_controls, arma_y_hat);

    // Compute Pearson correlation of the residuals
    partial_corr = arma::as_scalar(arma::cor(residuals_y, residuals_y_hat));

  } else {
    int i = controls.size();
    int j = controls.size() + 1;
    arma::mat data(clean_y.size(), i + 2);
    for (size_t i = 0; i < controls.size(); ++i) {
      data.col(i) = arma::vec(clean_controls[i]);
    }
    data.col(i) = arma::vec(clean_y);
    data.col(j) = arma::vec(clean_y_hat);

    // Compute the correlation matrix of the data
    arma::mat corrm = arma::cor(data);

    // Compute the precision matrix (inverse of the correlation matrix)
    arma::mat precm = arma::inv(corrm);

    // Get the correlation between y and y_hat after controlling for the others
    partial_corr = -precm(i, j) / std::sqrt(precm(i, i) * precm(j, j));
  }

  // Ensure partial correlation is within valid range [-1, 1]
  if (partial_corr < -1.0) partial_corr = -1.0;
  if (partial_corr > 1.0) partial_corr = 1.0;

  return partial_corr;
}

double PartialCorTrivar(const std::vector<double>& y,
                        const std::vector<double>& y_hat,
                        const std::vector<double>& control,
                        bool NA_rm = false,
                        bool linear = false){
  std::vector<std::vector<double>> conmat;
  conmat.push_back(control);

  double res = PartialCor(y,y_hat,conmat,NA_rm,linear);
  return res;
}

// Function to calculate the significance of a (partial) correlation coefficient
double CppCorSignificance(double r, int n, int k = 0) {
  double t = r * std::sqrt((n - k - 2) / (1 - r * r));
  return (1 - R::pt(t, n - 2, true, false)) * 2;
}

// Function to calculate the confidence interval for a (partial) correlation coefficient
std::vector<double> CppCorConfidence(double r, int n, int k = 0,
                                     double level = 0.05) {
  // Calculate the Fisher's z-transformation
  double z = 0.5 * std::log((1 + r) / (1 - r));

  // Calculate the standard error of z
  double ztheta = 1 / std::sqrt(n - k - 3);

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

// Function to find k-nearest neighbors of a given index in the embedding space
std::vector<std::size_t> CppKNNIndice(
    const std::vector<std::vector<double>>& embedding_space,
    std::size_t target_idx,
    std::size_t k)
{
  std::size_t n = embedding_space.size();
  std::vector<std::pair<double, std::size_t>> distances;

  for (std::size_t i = 0; i < n; ++i) {
    if (i == target_idx) continue;

    // Check if the entire embedding_space[i] is NaN
    if (std::all_of(embedding_space[i].begin(), embedding_space[i].end(),
                    [](double v) { return std::isnan(v); })) {
      continue;
    }

    double dist = CppDistance(embedding_space[target_idx], embedding_space[i], false, true);

    // Skip NaN distances
    if (!std::isnan(dist)) {
      distances.emplace_back(dist, i);
    }
  }

  // Partial sort to get k-nearest neighbors, excluding NaN distances
  std::partial_sort(distances.begin(), distances.begin() + std::min(k, distances.size()), distances.end());

  std::vector<std::size_t> neighbors;
  for (std::size_t i = 0; i < k && i < distances.size(); ++i) {
    neighbors.push_back(distances[i].second);
  }

  return neighbors;
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

// Function to remove linear trend using Armadillo for internal calculations
std::vector<double> LinearTrendRM(const std::vector<double>& vec,
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

// // Function to perform Linear Trend Removal
// std::vector<double> LinearTrendRM(const std::vector<double>& vec,
//                                   const std::vector<double>& xcoord,
//                                   const std::vector<double>& ycoord,
//                                   bool NA_rm = false) {
//   if (vec.size() != xcoord.size() || vec.size() != ycoord.size()) {
//     throw std::invalid_argument("Input vectors must have the same size.");
//   }
//
//   // Perform linear regression
//   double x1_mean = CppMean(xcoord, NA_rm);
//   double x2_mean = CppMean(ycoord, NA_rm);
//   double y_mean = CppMean(vec, NA_rm);
//
//   double x1_var = CppVariance(xcoord, NA_rm);
//   double x2_var = CppVariance(ycoord, NA_rm);
//   double x1_x2_cov = CppCovariance(xcoord, ycoord, NA_rm);
//
//   double x1_y_cov = CppCovariance(xcoord, vec, NA_rm);
//   double x2_y_cov = CppCovariance(ycoord, vec, NA_rm);
//
//   double denom = x1_var * x2_var - x1_x2_cov * x1_x2_cov;
//   if (denom == 0.0) {
//     throw std::invalid_argument("Linear regression cannot be performed due to collinearity.");
//   }
//
//   double b1 = (x2_var * x1_y_cov - x1_x2_cov * x2_y_cov) / denom;
//   double b2 = (x1_var * x2_y_cov - x1_x2_cov * x1_y_cov) / denom;
//   double b0 = y_mean - b1 * x1_mean - b2 * x2_mean;
//
//   // Predict vec_hat using the linear regression model
//   std::vector<double> vec_hat(vec.size());
//   for (size_t i = 0; i < vec.size(); ++i) {
//     vec_hat[i] = b0 + b1 * xcoord[i] + b2 * ycoord[i];
//   }
//
//   // Calculate vec - vec_hat
//   std::vector<double> result(vec.size());
//   for (size_t i = 0; i < vec.size(); ++i) {
//     result[i] = vec[i] - vec_hat[i];
//   }
//
//   return result;
// }

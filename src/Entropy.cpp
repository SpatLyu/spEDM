#include <cmath>
#include <vector>
#include <numeric>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "CppStats.h"

/**
 * @brief Computes the entropy of a given vector using k-nearest neighbors estimation.
 *
 * @param vec A vector of double values representing the dataset.
 * @param k The number of nearest neighbors to consider in the estimation.
 * @param base The logarithm base used for entropy calculation (default: 10).
 * @param NA_rm A boolean flag indicating whether to remove missing values (default: false).
 *
 * @return The estimated entropy of the vector.
 */
double CppEntropy_Cont(const std::vector<double>& vec, size_t k,
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

/**
 * @brief Computes the joint entropy of a multivariate matrix using k-nearest neighbors estimation.
 *
 * @param mat A 2D vector of double values where each row represents a data point.
 * @param columns The indices of columns to select for joint entropy calculation.
 * @param k The number of nearest neighbors to consider in the estimation.
 * @param base The logarithm base used for entropy calculation (default: 10).
 * @param NA_rm A boolean flag indicating whether to remove missing values (NaN) before computation (default: false).
 *
 * @return The estimated joint entropy of the multivariate matrix.
 */
double CppJoinEntropy_Cont(const std::vector<std::vector<double>>& mat,
                           const std::vector<int>& columns, size_t k,
                           double base = 10, bool NA_rm = false) {
  // Step 1: Construct new_mat based on selected columns
  std::vector<std::vector<double>> new_mat;
  size_t original_ncol = mat.empty() ? 0 : mat[0].size();

  // Check if columns match original matrix column count
  if (columns.size() == original_ncol) {
    new_mat = mat;  // Use original matrix directly if columns are identical
  } else {
    // Build new_mat by selecting specified columns
    new_mat.reserve(mat.size());
    for (const auto& row : mat) {
      std::vector<double> new_row;
      new_row.reserve(columns.size());
      for (int col : columns) {
        new_row.push_back(row[col]);
      }
      new_mat.push_back(std::move(new_row));
    }
  }

  // Step 2: Compute parameters based on new_mat
  size_t nrow = new_mat.size();
  size_t ncol = new_mat.empty() ? 0 : new_mat[0].size();

  // Step 3: Compute Chebyshev distance matrix for new_mat
  std::vector<double> distances(nrow);
  std::vector<std::vector<double>> mat_dist = CppMatChebyshevDistance(new_mat, NA_rm);

  // Step 4: Calculate k-th nearest neighbor distances
  for (size_t i = 0; i < nrow; ++i) {
    std::vector<double> dist_n;

    if (NA_rm) {
      for (double val : mat_dist[i]) {
        if (!std::isnan(val)) {
          dist_n.push_back(val);
        }
      }
    } else {
      dist_n = mat_dist[i];
    }

    // Handle k-th nearest neighbor selection
    if (k < dist_n.size()) {
      std::nth_element(dist_n.begin(), dist_n.begin() + k, dist_n.end());
      distances[i] = dist_n[k];
    } else {
      distances[i] = *std::max_element(dist_n.begin(), dist_n.end());
    }
  }

  // Step 5: Compute entropy components
  double sum = 0.0;
  for (size_t i = 0; i < nrow; i++) {
    sum += CppLog(2 * distances[i], base);
  }
  sum = sum * static_cast<double>(ncol) / nrow;

  // Final entropy calculation
  double E = CppDigamma(nrow) - CppDigamma(k) + sum;
  return E;
}

/**
 * @brief Computes the mutual information(MI) between two variables using k-nearest neighbors estimation.
 *
 * @note
 * True mutual information can't be negative. If its estimate by a numerical
 * method is negative, it means (providing the method is adequate) that the
 * mutual information is close to 0 and replacing it by 0 is a reasonable
 * strategy.
 *
 * @reference
 *  https://github.com/cran/NlinTS/blob/master/src/nsEntropy.cpp
 *  https://github.com/PengTao-HUST/crossmapy/blob/master/crossmapy/mi.py
 *
 * @param mat A 2D vector of double values where each row represents a data point with two variables.
 * @param k The number of nearest neighbors to consider in the estimation.
 * @param alg The algorithm choice for MI estimation (1: Kraskov Algorithm I, 2: Kraskov Algorithm II).
 * @param normalize A boolean flag indicating whether to normalize the MI by the joint entropy (default: false).
 * @param NA_rm A boolean flag indicating whether to remove missing values (NaN) before computation (default: false).
 *
 * @return The estimated mutual information between the two variables.
 */
double CppMutualInformation_Cont(const std::vector<std::vector<double>>& mat, size_t k,
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
  std::vector<std::vector<double>> mat_dist = CppMatChebyshevDistance(mat, NA_rm);

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

  // Mutual information is forced to 0 when it is negative
  mi = std::max(0.0, mi);

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

/**
 * @brief Computes the conditional entropy of x given y using k-nearest neighbors estimation.
 *
 * @param vecx A vector of double values representing the first variable.
 * @param vecy A vector of double values representing the second variable.
 * @param k The number of nearest neighbors to consider in the estimation.
 * @param base The logarithm base used for entropy calculation (default: 10).
 * @param NA_rm A boolean flag indicating whether to remove missing values (NaN) before computation (default: false).
 *
 * @return The estimated conditional entropy of x given y.
 */
double CppConditionalEntropy_Cont(const std::vector<double>& vecx,
                                  const std::vector<double>& vecy,
                                  size_t k, double base = 10,
                                  bool NA_rm = false) {
  // Create a 2D vector for the joint of x and y
  std::vector<std::vector<double>> joint_vec(vecx.size(), std::vector<double>(2));
  for (size_t i = 0; i < vecx.size(); ++i) {
    joint_vec[i][0] = vecx[i];
    joint_vec[i][1] = vecy[i];
  }

  // Compute the joint entropy H(X, Y)
  double joint_entropy = CppJoinEntropy_Cont(joint_vec, {0,1}, k, base, NA_rm);

  // Compute the entropy of y, H(Y)
  double entropy_y = CppEntropy_Cont(vecy, k, base, NA_rm);

  // Compute the conditional entropy H(X|Y) = H(X, Y) - H(Y)
  double ce = joint_entropy - entropy_y;

  return ce;
}

/**
 * Computes the entropy of a discrete sequence.
 * @param vec Input vector containing discrete values.
 * @param base Logarithm base (default: 10).
 * @param NA_rm If true, removes NaN values; otherwise returns NaN if any NaN exists.
 * @return Entropy value or NaN if invalid conditions occur.
 */
double CppEntropy_Disc(const std::vector<double>& vec,
                       double base = 10, bool NA_rm = false) {
  std::unordered_map<double, int> counts;
  int valid_n = 0;

  for (double x : vec) {
    if (std::isnan(x)) {
      if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    counts[x]++;
    valid_n++;
  }

  if (valid_n == 0) return std::numeric_limits<double>::quiet_NaN();

  const double log_base = std::log(base);
  double entropy = 0.0;

  for (const auto& pair : counts) {
    double p = static_cast<double>(pair.second) / valid_n;
    entropy += p * std::log(p);
  }

  return -entropy / log_base;
}

/**
 * Computes the joint entropy of a multivariate discrete sequence.
 * @param mat Input matrix where each row represents a sample containing multiple variables.
 * @param columns The columns which used in joint entropy estimation.
 * @param base Logarithm base (default: 10).
 * @param NA_rm If true, removes samples with any NaN; otherwise returns NaN if any NaN exists.
 * @return Joint entropy value or NaN if invalid conditions occur.
 */
double CppJoinEntropy_Disc(const std::vector<std::vector<double>>& mat,
                           const std::vector<int>& columns,
                           double base = 10, bool NA_rm = false){
  const double log_base = std::log(base);

  // Flattened and valid samples, stored as string key or unique encoding
  std::unordered_map<std::string, int> counts;
  int valid_count = 0;

  for (const auto& sample : mat) {
    bool has_nan = false;
    std::string key;

    for (size_t i = 0; i < columns.size(); ++i) {
      double val = sample[columns[i]];
      if (std::isnan(val)) {
        has_nan = true;
        break;
      }
      // simple separator-based encoding
      key += std::to_string(val) + "_";
    }

    if (has_nan) {
      if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    counts[key]++;
    valid_count++;
  }

  if (valid_count == 0) return std::numeric_limits<double>::quiet_NaN();

  double entropy = 0.0;
  for (const auto& pair : counts) {
    double p = static_cast<double>(pair.second) / valid_count;
    entropy += p * std::log(p);
  }

  return -entropy / log_base;
}

/**
 * Computes the mutual information between two sets of discrete variables (columns).
 * @param mat Input matrix where each row represents a sample and each column a discrete variable.
 * @param columns1 Indices of columns representing the first set of variables (X).
 * @param columns2 Indices of columns representing the second set of variables (Y).
 * @param base Logarithm base used in entropy calculations (default: 10).
 * @param NA_rm If true, removes samples with any NaN values; otherwise returns NaN if any NaN is encountered.
 * @return Mutual information value I(X; Y) = H(X) + H(Y) - H(X,Y), or NaN if invalid conditions occur.
 */
double CppMutualInformation_Disc(const std::vector<std::vector<double>>& mat,
                                 const std::vector<int>& columns1,
                                 const std::vector<int>& columns2,
                                 double base = 10, bool NA_rm = false) {
  std::unordered_set<int> unique_set;
  unique_set.insert(columns1.begin(), columns1.end());
  unique_set.insert(columns2.begin(), columns2.end());
  std::vector<int> columns(unique_set.begin(), unique_set.end());

  double h_x = CppJoinEntropy_Disc(mat, columns1, base, NA_rm);
  double h_y = CppJoinEntropy_Disc(mat, columns2, base, NA_rm);
  double h_xy = CppJoinEntropy_Disc(mat, columns, base, NA_rm);

  if (std::isnan(h_x) || std::isnan(h_y) || std::isnan(h_xy)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return h_x + h_y - h_xy;
}

/**
 * Computes the conditional entropy H(X | Y) between two sets of discrete variables.
 * @param mat Input matrix where each row is a sample and each column is a discrete variable.
 * @param target_columns Indices of columns representing the target variable(s) X.
 * @param conditional_columns Indices of columns representing the conditioning variable(s) Y.
 * @param base Logarithm base used in entropy calculations (default: 10).
 * @param NA_rm If true, removes samples with any NaN values; otherwise returns NaN if any NaN is encountered.
 * @return Conditional entropy value H(X | Y) = H(X,Y) - H(Y), or NaN if invalid conditions occur.
 */
double CppConditionalEntropy_Disc(const std::vector<std::vector<double>>& mat,
                                  const std::vector<int>& target_columns,
                                  const std::vector<int>& conditional_columns,
                                  double base = 10, bool NA_rm = false) {
  std::unordered_set<int> unique_set;
  unique_set.insert(target_columns.begin(), target_columns.end());
  unique_set.insert(conditional_columns.begin(), conditional_columns.end());
  std::vector<int> columns(unique_set.begin(), unique_set.end());

  double H_xy = CppJoinEntropy_Disc(mat, columns, base, NA_rm);
  double H_y = CppJoinEntropy_Disc(mat, conditional_columns, base, NA_rm);

  if (std::isnan(H_xy) || std::isnan(H_y)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return H_xy - H_y; // H(X|Y) = H(X,Y) - H(Y)
}

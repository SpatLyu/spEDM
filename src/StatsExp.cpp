#include <vector>
#include "CppStats.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double RcppMean(const Rcpp::NumericVector& vec,
                bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppMean(y, NA_rm);
}

// [[Rcpp::export]]
double RcppSum(const Rcpp::NumericVector& vec,
               bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppSum(y, NA_rm);
}

// [[Rcpp::export]]
double RcppVariance(const Rcpp::NumericVector& vec,
                    bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppVariance(y, NA_rm);
}

// [[Rcpp::export]]
double RcppCovariance(const Rcpp::NumericVector& vec1,
                      const Rcpp::NumericVector& vec2,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppMAE function
  return CppCovariance(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
double RcppMAE(const Rcpp::NumericVector& vec1,
               const Rcpp::NumericVector& vec2,
               bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppMAE function
  return CppMAE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
double RcppRMSE(const Rcpp::NumericVector& vec1,
                const Rcpp::NumericVector& vec2,
                bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppRMSE function
  return CppRMSE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
Rcpp::NumericVector RcppAbsDiff(const Rcpp::NumericVector& vec1,
                                const Rcpp::NumericVector& vec2) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec1_cpp = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> vec2_cpp = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppAbsDiff function
  std::vector<double> result = CppAbsDiff(vec1_cpp, vec2_cpp);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
Rcpp::NumericVector RcppSumNormalize(const Rcpp::NumericVector& vec, bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_cpp = Rcpp::as<std::vector<double>>(vec);

  // Call the CppSumNormalize function
  std::vector<double> result = CppSumNormalize(vec_cpp, NA_rm);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
double RcppDistance(const Rcpp::NumericVector& vec1,
                    const Rcpp::NumericVector& vec2,
                    bool L1norm = false,
                    bool NA_rm = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> v2 = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppDistance function
  return CppDistance(v1, v2, L1norm ,NA_rm);
}

// [[Rcpp::export]]
double RcppPearsonCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_vec = Rcpp::as<std::vector<double>>(y);
  std::vector<double> y_hat_vec = Rcpp::as<std::vector<double>>(y_hat);

  // Call the PearsonCor function
  return PearsonCor(y_vec, y_hat_vec, NA_rm);
}

// Rcpp wrapper for PartialCor function
// [[Rcpp::export]]
double RcppPartialCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      const Rcpp::NumericMatrix& controls,
                      bool NA_rm = false,
                      bool linear = false) {

  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> std_controls(controls.ncol());
  for (int i = 0; i < controls.ncol(); ++i) {
    Rcpp::NumericVector covvar = controls.column(i);
    std_controls[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Call the PartialCor function
  return PartialCor(std_y, std_y_hat, std_controls, NA_rm, linear);
}

// [[Rcpp::export]]
double RcppPartialCorTrivar(const Rcpp::NumericVector& y,
                            const Rcpp::NumericVector& y_hat,
                            const Rcpp::NumericVector& control,
                            bool NA_rm = false,
                            bool linear = false) {

  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);
  std::vector<double> std_control = Rcpp::as<std::vector<double>>(control);

  // Call the PartialCorTrivar function
  return PartialCorTrivar(std_y, std_y_hat, std_control, NA_rm, linear);
}

// Wrapper function to calculate the significance of a (partial) correlation coefficient
// [[Rcpp::export]]
double RcppCorSignificance(double r, int n, int k = 0){
  return CppCorSignificance(r, n, k);
}

// Wrapper function to calculate the confidence interval for a (partial) correlation coefficient and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppCorConfidence(double r, int n, int k = 0, double level = 0.05) {
  // Calculate the confidence interval
  std::vector<double> result = CppCorConfidence(r, n, k, level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to compute the distance matrix of a given matrix 'mat'
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppMatDistance(const Rcpp::NumericMatrix& mat,
                                    bool L1norm = false, bool NA_rm = false) {

  // Convert the Rcpp::NumericMatrix to a C++ vector of vectors (std::vector)
  size_t n = mat.nrow();
  std::vector<std::vector<double>> cppMat(n, std::vector<double>(mat.ncol()));

  // Fill cppMat with values from the R matrix
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < mat.ncol(); ++j) {
      cppMat[i][j] = mat(i, j);
    }
  }

  // Call the C++ function to compute the distance matrix
  std::vector<std::vector<double>> distanceMatrix = CppMatDistance(cppMat, L1norm, NA_rm);

  // Convert the resulting distance matrix back into an Rcpp::NumericMatrix
  Rcpp::NumericMatrix result(n, n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      result(i, j) = distanceMatrix[i][j];
    }
  }

  return result;
}

// Wrapper function to find k-nearest neighbors of a given index in the embedding space
// [[Rcpp::export]]
Rcpp::IntegerVector RcppKNNIndice(const Rcpp::NumericMatrix& embedding_space,
                                  int target_idx, int k) {
  // Get the number of rows and columns
  std::size_t n_rows = embedding_space.nrow();
  std::size_t n_cols = embedding_space.ncol();

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> embedding_vec(n_rows, std::vector<double>(n_cols));
  for (std::size_t i = 0; i < n_rows; ++i) {
    for (std::size_t j = 0; j < n_cols; ++j) {
      embedding_vec[i][j] = embedding_space(i, j);
    }
  }

  // Ensure target index is within valid range
  if (target_idx < 0 || static_cast<std::size_t>(target_idx) >= n_rows) {
    Rcpp::stop("target_idx is out of range.");
  }

  // Ensure k is positive
  if (k <= 0) {
    Rcpp::stop("k must be greater than 0.");
  }

  // Call the C++ function
  std::vector<std::size_t> knn_indices = CppKNNIndice(embedding_vec,
                                                      static_cast<std::size_t>(target_idx - 1),
                                                      static_cast<std::size_t>(k));

  // Convert result to Rcpp::IntegerVector (R uses 1-based indexing)
  Rcpp::IntegerVector result(knn_indices.size());
  for (std::size_t i = 0; i < knn_indices.size(); ++i) {
    result[i] = static_cast<int>(knn_indices[i]) + 1;  // Convert to 1-based index
  }

  return result;
}

// Wrapper function to find k-nearest neighbors of a given index using a precomputed distance matrix
// [[Rcpp::export]]
Rcpp::IntegerVector RcppDistKNNIndice(const Rcpp::NumericMatrix& dist_mat,
                                      int target_idx, int k) {
  // Get the number of rows and columns
  std::size_t n_rows = dist_mat.nrow();
  std::size_t n_cols = dist_mat.ncol();

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> distmat(n_rows, std::vector<double>(n_cols));
  for (std::size_t i = 0; i < n_rows; ++i) {
    for (std::size_t j = 0; j < n_cols; ++j) {
      distmat[i][j] = dist_mat(i, j);
    }
  }

  // Ensure target index is within valid range
  if (target_idx < 0 || static_cast<std::size_t>(target_idx) >= n_rows) {
    Rcpp::stop("target_idx is out of range.");
  }

  // Ensure k is positive
  if (k <= 0) {
    Rcpp::stop("k must be greater than 0.");
  }

  // Call the C++ function
  std::vector<std::size_t> knn_indices = CppDistKNNIndice(distmat,
                                                          static_cast<std::size_t>(target_idx - 1),
                                                          static_cast<std::size_t>(k));

  // Convert result to Rcpp::IntegerVector (R uses 1-based indexing)
  Rcpp::IntegerVector result(knn_indices.size());
  for (std::size_t i = 0; i < knn_indices.size(); ++i) {
    result[i] = static_cast<int>(knn_indices[i]) + 1;  // Convert to 1-based index
  }

  return result;
}

// Wrapper function to perform Linear Trend Removal and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppLinearTrendRM(const Rcpp::NumericVector& vec,
                                      const Rcpp::NumericVector& xcoord,
                                      const Rcpp::NumericVector& ycoord,
                                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
  std::vector<double> xcoord_std = Rcpp::as<std::vector<double>>(xcoord);
  std::vector<double> ycoord_std = Rcpp::as<std::vector<double>>(ycoord);

  // Perform Linear Trend Removal
  std::vector<double> result = LinearTrendRM(vec_std, xcoord_std, ycoord_std, NA_rm);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Rcpp wrapper function for CppSVD
// [[Rcpp::export]]
Rcpp::List RcppSVD(const Rcpp::NumericMatrix& X) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t m = X.nrow();
  size_t n = X.ncol();
  std::vector<std::vector<double>> X_vec(m, std::vector<double>(n));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      X_vec[i][j] = X(i, j);
    }
  }

  // Call the original CppSVD function
  std::vector<std::vector<std::vector<double>>> result = CppSVD(X_vec);

  // Extract results from CppSVD output
  std::vector<std::vector<double>> u = result[0]; // Left singular vectors
  std::vector<double> d = result[1][0];           // Singular values
  std::vector<std::vector<double>> v = result[2]; // Right singular vectors

  // Convert std::vector results to Rcpp objects
  Rcpp::NumericMatrix u_rcpp(m, m);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      u_rcpp(i, j) = u[i][j];
    }
  }

  Rcpp::NumericVector d_rcpp(d.size());
  for (size_t i = 0; i < d.size(); ++i) {
    d_rcpp(i) = d[i];
  }

  Rcpp::NumericMatrix v_rcpp(v.size(), v[0].size());
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v_rcpp(i, j) = v[i][j];
    }
  }

  // Return results as an Rcpp::List to match R's svd() output
  return Rcpp::List::create(
    Rcpp::Named("u") = u_rcpp, // Left singular vectors
    Rcpp::Named("d") = d_rcpp, // Singular values
    Rcpp::Named("v") = v_rcpp  // Right singular vectors
  );
}

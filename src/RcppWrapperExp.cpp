#include <Rcpp.h>
#include "CppStats.h"
#include "CppUtils.h"

// Wrapper function to calculate the confidence interval for a correlation coefficient and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppConfidence(double r, int n, double level = 0.05) {
  // Calculate the confidence interval
  std::vector<double> result = CppConfidence(r, n, level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
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

// Wrapper function to calculate lagged indices and return a List
// [[Rcpp::export]]
Rcpp::List RcppLaggedIndices(const Rcpp::NumericVector& vec,
                             const Rcpp::NumericMatrix& nbmat,
                             int lagNum) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Calculate lagged indices
  std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec_std, nbmat_std, lagNum);

  // Convert std::vector<std::vector<int>> to Rcpp::List
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = Rcpp::wrap(lagged_indices[i]);
  }

  return result;
}

// Wrapper function to generate embeddings and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGenEmbeddings(const Rcpp::NumericVector& vec,
                                      const Rcpp::NumericMatrix& nbmat,
                                      int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenEmbeddings(vec_std, nbmat_std, E);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  int rows = embeddings.size();
  int cols = embeddings[0].size();
  Rcpp::NumericMatrix result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = embeddings[i][j];
    }
  }

  return result;
}

#include <Rcpp.h>
#include "CppStats.h"
#include "CppUtils.h"
#include "CppGCCM.h"

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

// Wrapper function to calculate the pairwise absolute difference mean matrix
// and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppDist(const Rcpp::NumericMatrix& matrix) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t n = matrix.nrow();
  size_t m = matrix.ncol();
  std::vector<std::vector<double>> matrix_std(n, std::vector<double>(m));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      matrix_std[i][j] = matrix(i, j);
    }
  }

  // Calculate the pairwise absolute difference mean matrix
  std::vector<std::vector<double>> dist_matrix = CppDist(matrix_std);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix result(n, n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      result(i, j) = dist_matrix[i][j];
    }
  }

  return result;
}

// Wrapper function to find the indices of the libsize + 1 closest elements for each row in distmat
// and return an IntegerMatrix
// [[Rcpp::export]]
Rcpp::IntegerMatrix RcppClosestIndices(const Rcpp::NumericMatrix& distmat, int libsize) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t n = distmat.nrow();
  std::vector<std::vector<double>> distmat_std(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      distmat_std[i][j] = distmat(i, j);
    }
  }

  // Find the indices of the libsize + 1 closest elements for each row in distmat
  std::vector<std::vector<int>> closestIndices = CppClosestIndices(distmat_std, libsize);

  // Convert std::vector<std::vector<int>> to Rcpp::IntegerMatrix
  Rcpp::IntegerMatrix result(n, libsize + 1);
  for (size_t i = 0; i < n; ++i) {
    for (int j = 0; j < libsize + 1; ++j) {
      result(i, j) = closestIndices[i][j];
    }
  }

  return result;
}

// Wrapper function to calculate the CCM weights and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppCCMWeight(const Rcpp::NumericMatrix& distmat,
                                  const Rcpp::IntegerMatrix& closestIndices,
                                  int libsize) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t n = distmat.nrow();
  std::vector<std::vector<double>> distmat_std(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      distmat_std[i][j] = distmat(i, j);
    }
  }

  // Convert Rcpp::IntegerMatrix to std::vector<std::vector<int>>
  std::vector<std::vector<int>> closestIndices_std(n, std::vector<int>(libsize + 1));
  for (size_t i = 0; i < n; ++i) {
    for (int j = 0; j < libsize + 1; ++j) {
      closestIndices_std[i][j] = closestIndices(i, j);
    }
  }

  // Calculate the CCM weights
  std::vector<std::vector<double>> weights = CCMWeight(distmat_std, closestIndices_std, libsize);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix result(n, libsize + 1);
  for (size_t i = 0; i < n; ++i) {
    for (int j = 0; j < libsize + 1; ++j) {
      result(i, j) = weights[i][j];
    }
  }

  return result;
}

// Wrapper function to perform Simplex Projection and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppSimplexProjection(const Rcpp::NumericVector& y,
                                          const Rcpp::NumericVector& x,
                                          const Rcpp::NumericMatrix& nbmat,
                                          int libsize,
                                          int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Perform Simplex Projection
  std::vector<double> result = SimplexProjection(y_std, x_std, nbmat_std, libsize, E);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to perform GCCM Lattice and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCCMLattice(const Rcpp::NumericVector& y,
                                    const Rcpp::NumericVector& x,
                                    const Rcpp::NumericMatrix& nbmat,
                                    const Rcpp::IntegerVector& libsizes,
                                    int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);

  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = GCCMLattice(y_std, x_std, nbmat_std, libsizes_std, E);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 5);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
  }

  return resultMatrix;
}

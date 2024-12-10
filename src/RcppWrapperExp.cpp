#include <Rcpp.h>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "GCCM4Lattice.h"

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
Rcpp::List RcppLaggedVar4Lattice(const Rcpp::List& nb, int lagNum) {
  int n = nb.size();

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Calculate lagged indices
  std::vector<std::vector<int>> lagged_indices = CppLaggedVar4Lattice(nb_vec, lagNum);

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
                                      const Rcpp::List& nb,
                                      int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenEmbeddings(vec_std, nb_vec, E);

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

// Wrapper function to perform GCCM Lattice and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCCMLattice(const Rcpp::NumericVector& x,
                                    const Rcpp::NumericVector& y,
                                    const Rcpp::List& nb,
                                    const Rcpp::IntegerVector& libsizes,
                                    int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenEmbeddings(x_std, nb_vec, E);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);

  // Define the interval [1, n] as a std::vector<std::pair<int, int>>
  int n = nb_vec.size();
  std::vector<std::pair<int, int>> interval = {{1, n}};

  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = GCCMLattice(embeddings, y_std, libsizes_std,
                                                        interval, interval, E);

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

// // Wrapper function to calculate lagged indices and return a List
// // [[Rcpp::export]]
// Rcpp::List RcppLaggedIndices(const Rcpp::NumericVector& vec,
//                              const Rcpp::NumericMatrix& nbmat,
//                              int lagNum) {
//   // Convert Rcpp::NumericVector to std::vector<double>
//   std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
//
//   // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
//   int n = nbmat.nrow();
//   std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < n; ++j) {
//       nbmat_std[i][j] = nbmat(i, j);
//     }
//   }
//
//   // Calculate lagged indices
//   std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec_std, nbmat_std, lagNum);
//
//   // Convert std::vector<std::vector<int>> to Rcpp::List
//   Rcpp::List result(n);
//   for (int i = 0; i < n; ++i) {
//     result[i] = Rcpp::wrap(lagged_indices[i]);
//   }
//
//   return result;
// }
//
// // Wrapper function to generate embeddings and return a NumericMatrix
// // [[Rcpp::export]]
// Rcpp::NumericMatrix RcppGenEmbeddings(const Rcpp::NumericVector& vec,
//                                       const Rcpp::NumericMatrix& nbmat,
//                                       int E) {
//   // Convert Rcpp::NumericVector to std::vector<double>
//   std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
//
//   // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
//   int n = nbmat.nrow();
//   std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < n; ++j) {
//       nbmat_std[i][j] = nbmat(i, j);
//     }
//   }
//
//   // Generate embeddings
//   std::vector<std::vector<double>> embeddings = GenEmbeddings(vec_std, nbmat_std, E);
//
//   // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
//   int rows = embeddings.size();
//   int cols = embeddings[0].size();
//   Rcpp::NumericMatrix result(rows, cols);
//   for (int i = 0; i < rows; ++i) {
//     for (int j = 0; j < cols; ++j) {
//       result(i, j) = embeddings[i][j];
//     }
//   }
//
//   return result;
// }

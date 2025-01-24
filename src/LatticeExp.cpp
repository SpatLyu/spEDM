#include <vector>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "GCCM4Lattice.h"
#include <Rcpp.h>

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

// Rcpp wrapper function for ArmaLinearTrendRM
// [[Rcpp::export]]
Rcpp::NumericVector RcppArmaLinearTrendRM(const Rcpp::NumericVector& vec,
                                          const Rcpp::NumericVector& xcoord,
                                          const Rcpp::NumericVector& ycoord,
                                          bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_cpp(vec.begin(), vec.end());
  std::vector<double> xcoord_cpp(xcoord.begin(), xcoord.end());
  std::vector<double> ycoord_cpp(ycoord.begin(), ycoord.end());

  // Call the original ArmaLinearTrendRM function
  std::vector<double> result = ArmaLinearTrendRM(vec_cpp, xcoord_cpp, ycoord_cpp, NA_rm);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::NumericVector(result.begin(), result.end());
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

// Wrapper function to generate embeddings and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGenLatticeEmbeddings(const Rcpp::NumericVector& vec,
                                             const Rcpp::List& nb,
                                             int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec_std, nb_vec, E);

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
Rcpp::NumericMatrix RcppGCCM4Lattice(const Rcpp::NumericVector& x,
                                     const Rcpp::NumericVector& y,
                                     const Rcpp::List& nb,
                                     const Rcpp::IntegerVector& libsizes,
                                     int E,
                                     int tau,
                                     int b,
                                     bool simplex,
                                     double theta,
                                     bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(x_std, nb_vec, E);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);

  // Define the interval [1, n] as a std::vector<std::pair<int, int>>
  int n = nb_vec.size();
  std::vector<std::pair<int, int>> interval = {{1, n}};

  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = GCCM4Lattice(
    embeddings,
    y_std,
    libsizes_std,
    interval,
    interval,
    E,
    tau,
    b,
    simplex,
    theta,
    progressbar);

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

#include <vector>
#include "CppLatticeUtils.h"
#include "Forecast4Lattice.h"
#include "GCCM4Lattice.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>

// Function to convert Rcpp::List to std::vector<std::vector<int>>
std::vector<std::vector<int>> nb2vec(Rcpp::List nb) {
  // Get the number of elements in the nb object
  int n = nb.size();

  // Create a vector<vector<int>> to store the result
  std::vector<std::vector<int>> result(n);

  // Iterate over each element in the nb object
  for (int i = 0; i < n; ++i) {
    // Get the current element (should be an integer vector)
    Rcpp::IntegerVector current_nb = nb[i];

    // Create a vector<int> to store the current subset of elements
    std::vector<int> current_subset;

    // Iterate over each element in the current subset
    for (int j = 0; j < current_nb.size(); ++j) {
      // Subtract one from each element to convert from R's 1-based indexing to C++'s 0-based indexing
      current_subset.push_back(current_nb[j] - 1);
    }

    // Add the current subset to the result
    result[i] = current_subset;
  }

  return result;
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
Rcpp::NumericMatrix RcppGenLatticeEmbeddings(const Rcpp::NumericVector& vec,
                                             const Rcpp::List& nb,
                                             int E,
                                             bool includeself) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec_std, nb_vec, E, includeself);

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

// Description: Computes Simplex projection for lattice data and returns a matrix containing
//              the embedding dimension (E), Pearson correlation coefficient (PearsonCor),
//              mean absolute error (MAE), and root mean squared error (RMSE).
// Parameters:
//   - x: A NumericVector containing the time series data.
//   - nb: A List containing neighborhood information for lattice data.
//   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
//   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
//   - E: An IntegerVector specifying the embedding dimensions to test.
//   - b: An integer specifying the number of neighbors to use for simplex projection.
// Returns: A NumericMatrix where each row contains {E, PearsonCor, CppMAE, CppRMSE}.
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSimplex4Lattice(const Rcpp::NumericVector& x,
                                        const Rcpp::List& nb,
                                        const Rcpp::IntegerVector& lib,
                                        const Rcpp::IntegerVector& pred,
                                        const Rcpp::IntegerVector& E,
                                        int b,
                                        int threads) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(x);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(vec_std.size(), false);
  std::vector<bool> pred_indices(vec_std.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  for (int i = 0; i < lib.size(); ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
  }

  std::vector<std::vector<double>> res_std = Simplex4Lattice(vec_std,nb_vec,lib_indices,pred_indices,E_std,b,threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "rho", "mae", "rmse");
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
                                     int threads,
                                     bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(x_std, nb_vec, E, false);

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
    threads,
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

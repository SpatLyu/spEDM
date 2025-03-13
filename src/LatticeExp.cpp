#include <vector>
#include <cmath>
#include <algorithm>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "Forecast4Lattice.h"
#include "MultiViewEmbedding.h"
#include "GCCM4Lattice.h"
#include "SCPCM4Lattice.h"
#include "CrossMappingCardinality.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>

// Function to convert Rcpp::List to std::vector<std::vector<int>>
std::vector<std::vector<int>> nb2vec(const Rcpp::List& nb) {
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

// Wrapper function to calculate accumulated lagged neighbor indices and return a List
// [[Rcpp::export]]
Rcpp::List RcppLaggedNeighbor4Lattice(const Rcpp::List& nb, int lagNum) {
  int n = nb.size();

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Calculate lagged indices
  std::vector<std::vector<int>> lagged_indices = CppLaggedNeighbor4Lattice(nb_vec, lagNum);
  // Restore the 0-based index from C++ to the 1-based index in R
  for (auto& row : lagged_indices) {
    for (auto& val : row) {
      val += 1;
    }
  }
  // // A more modern C++ implementation, provided for comparison only. -- Wenbo Lv, 2025.3.3
  // for (auto& row : lagged_indices) {
  //   std::transform(row.begin(), row.end(), row.begin(), [](int x) { return x + 1; });
  // }

  // Convert std::vector<std::vector<int>> to Rcpp::List
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = Rcpp::wrap(lagged_indices[i]);
  }

  return result;
}

// Wrapper function to calculate lagged values and return a List
// [[Rcpp::export]]
Rcpp::List RcppLaggedVar4Lattice(const Rcpp::NumericVector& vec,
                                 const Rcpp::List& nb, int lagNum) {
  int n = nb.size();

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Calculate lagged indices
  std::vector<std::vector<double>> lagged_values = CppLaggedVar4Lattice(vec_std, nb_vec, lagNum);

  // Convert std::vector<std::vector<int>> to Rcpp::List
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = Rcpp::wrap(lagged_values[i]);
  }

  return result;
}

// Wrapper function to generate embeddings and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGenLatticeEmbeddings(const Rcpp::NumericVector& vec,
                                             const Rcpp::List& nb,
                                             int E,
                                             int tau) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec_std, nb_vec, E, tau);

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

/**
 * Description:
 *   Computes Simplex projection for lattice data and returns a matrix containing
 *   the embedding dimension (E), Pearson correlation coefficient (PearsonCor),
 *   mean absolute error (MAE), and root mean squared error (RMSE).
 *
 * Parameters:
 *   - x: A NumericVector containing the time series data.
 *   - nb: A List containing neighborhood information for lattice data.
 *   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
 *   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 *   - E: An IntegerVector specifying the embedding dimensions to test.
 *   - b: An IntegerVector specifying the numbers of neighbors to use for simplex projection.
 *   - tau: An integer specifying the step of spatial lags for prediction.
 *
 * Returns:
 *   A NumericMatrix where each row contains {E, b, PearsonCor, MAE, RMSE}:
 *   - E: The tested embedding dimension.
 *   - b: The tested numbers of neighbors
 *   - PearsonCor: The Pearson correlation coefficient between the predicted and actual values.
 *   - MAE: The mean absolute error between the predicted and actual values.
 *   - RMSE: The root mean squared error between the predicted and actual values.
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSimplex4Lattice(const Rcpp::NumericVector& x,
                                        const Rcpp::List& nb,
                                        const Rcpp::IntegerVector& lib,
                                        const Rcpp::IntegerVector& pred,
                                        const Rcpp::IntegerVector& E,
                                        const Rcpp::IntegerVector& b,
                                        int tau,
                                        int threads) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(x);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(vec_std.size(), false);
  std::vector<bool> pred_indices(vec_std.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
  }

  std::vector<std::vector<double>> res_std = Simplex4Lattice(
    vec_std,
    nb_vec,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    threads);

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
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "rho", "mae", "rmse");
  return result;
}

/**
 * Parameters:
 *   - x: A NumericVector containing the time series data.
 *   - nb: A List containing neighborhood information for lattice data.
 *   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
 *   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 *   - theta: A NumericVector containing the parameter values to be tested for theta.
 *   - E: An integer specifying the embedding dimension to test.
 *   - tau: the step of spatial lags for prediction.
 *   - b: An integer specifying the number of neighbors to use for s-mapping.
 *   - threads: An integer specifying the number of threads to use for parallel computation.
 *
 * Returns:
 *   A NumericMatrix where each row contains {theta, PearsonCor, MAE, RMSE}:
 *   - theta: The tested parameter value.
 *   - PearsonCor: The Pearson correlation coefficient between the predicted and actual values.
 *   - MAE: The mean absolute error between the predicted and actual values.
 *   - RMSE: The root mean squared error between the predicted and actual values.
 *
 * This function utilizes parallel processing for faster computation across different values of theta.
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSMap4Lattice(const Rcpp::NumericVector& x,
                                     const Rcpp::List& nb,
                                     const Rcpp::IntegerVector& lib,
                                     const Rcpp::IntegerVector& pred,
                                     const Rcpp::NumericVector& theta,
                                     int E,
                                     int tau,
                                     int b,
                                     int threads) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(vec_std.size(), false);
  std::vector<bool> pred_indices(vec_std.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
  }

  std::vector<std::vector<double>> res_std = SMap4Lattice(
    vec_std,
    nb_vec,
    lib_indices,
    pred_indices,
    theta_std,
    E,
    tau,
    b,
    threads);

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
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("theta", "rho", "mae", "rmse");
  return result;
}

/*
 * Wrapper function to perform multiview embedding for spatial lattice data.
 *
 * Parameters:
 * - x: An Rcpp::NumericMatrix representing the selected multivariate spatial cross sectional data input.
 * - y: An Rcpp::NumericVector representing the target variable for prediction.
 * - nb: An Rcpp::List representing the spatial neighborhood structure.
 * - lib: An Rcpp::IntegerVector specifying the library indices (1-based, from R).
 * - pred: An Rcpp::IntegerVector specifying the prediction indices (1-based, from R).
 * - E: An integer specifying the embedding dimensions.
 * - tau: An integer specifying the step of spatial lags.
 * - b: An integer specifying the number of neighbors to use for simplex projection.
 * - top: An integer specifying the number of top embeddings to consider; if <= 0, uses sqrt(m) heuristic.
 * - nvar: An integer specifying the number of `nvar`-dimensional variable combinations.
 * - threads: An integer indicating the number of threads for parallel processing.
 *
 * Returns:
 * - An Rcpp::NumericVector containing the prediction results based on the multiview embedding.
 */
// [[Rcpp::export]]
Rcpp::NumericVector RcppMultiView4Lattice(const Rcpp::NumericMatrix& x,
                                          const Rcpp::NumericVector& y,
                                          const Rcpp::List& nb,
                                          const Rcpp::IntegerVector& lib,
                                          const Rcpp::IntegerVector& pred,
                                          int E,
                                          int tau,
                                          int b,
                                          int top,
                                          int nvar,
                                          int threads){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> target = Rcpp::as<std::vector<double>>(y);

  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(target.size(), false);
  std::vector<bool> pred_indices(target.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
  }

  int num_row = x.nrow();
  int num_var = x.ncol();

  //  if top <= 0, we choose to apply the heuristic of k (sqrt(m))
  int k;
  if (top <= 0){
    double m = CppCombine(num_var*E,nvar) - CppCombine(num_var*(E - 1),nvar);
    k = std::floor(std::sqrt(m));
  } else {
    k = top;
  }

  // Combine all the lags in the embeddings
  std::vector<std::vector<double>> vec_std(num_row,std::vector<double>(E*num_var,std::numeric_limits<double>::quiet_NaN()));
  for (int n = 0; n < num_var; ++n) {
    // Initialize a std::vector to store the column values
    std::vector<double> univec(num_row);

    // Copy the nth column from the matrix to the vector
    for (int i = 0; i < num_row; ++i) {
      univec[i] = x(i, n);  // Access element at (i, n)
    }

    // Generate the embedding:
    std::vector<std::vector<double>> vectors = GenLatticeEmbeddings(univec,nb_vec,E,tau);

    for (size_t row = 0; row < vectors.size(); ++row) {  // Loop through each row
      for (size_t col = 0; col < vectors[0].size(); ++col) {  // Loop through each column
        vec_std[row][n * E + col] = vectors[row][col];  // Copy elements
      }
    }
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns; // To store indices of valid columns

  // Iterate over each column to check if it contains any non-NaN values
  for (size_t col = 0; col < vec_std[0].size(); ++col) {
    bool isAllNaN = true;
    for (size_t row = 0; row < vec_std.size(); ++row) {
      if (!std::isnan(vec_std[row][col])) {
        isAllNaN = false;
        break;
      }
    }
    if (!isAllNaN) {
      validColumns.push_back(col); // Store the index of valid columns
    }
  }

  if (validColumns.size() != vec_std[0].size()) {
    std::vector<std::vector<double>> filteredEmbeddings;
    for (size_t row = 0; row < vec_std.size(); ++row) {
      std::vector<double> filteredRow;
      for (size_t col : validColumns) {
        filteredRow.push_back(vec_std[row][col]);
      }
      filteredEmbeddings.push_back(filteredRow);
    }
    vec_std = filteredEmbeddings;
  }

  std::vector<double> res = MultiViewEmbedding(
    vec_std,
    target,
    lib_indices,
    pred_indices,
    b,
    k,
    threads);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(res);
}

// Wrapper function to perform GCCM Lattice and return a NumericMatrix
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCCM4Lattice(const Rcpp::NumericVector& x,
                                     const Rcpp::NumericVector& y,
                                     const Rcpp::List& nb,
                                     const Rcpp::IntegerVector& libsizes,
                                     const Rcpp::IntegerVector& lib,
                                     const Rcpp::IntegerVector& pred,
                                     int E,
                                     int tau,
                                     int b,
                                     bool simplex,
                                     double theta,
                                     int threads,
                                     int parallel_level,
                                     bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);
  std::vector<int> lib_std = Rcpp::as<std::vector<int>>(lib);
  std::vector<int> pred_std = Rcpp::as<std::vector<int>>(pred);

  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = GCCM4Lattice(
    x_std,
    y_std,
    nb_vec,
    libsizes_std,
    lib_std,
    pred_std,
    E,
    tau,
    b,
    simplex,
    theta,
    threads,
    parallel_level,
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

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("libsizes",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_upper","x_xmap_y_lower");
  return resultMatrix;
}

// Wrapper function to perform SCPCM Lattice and return a NumericMatrix
// predict y based on x ====> x xmap y ====> y causes x (account for controls)
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSCPCM4Lattice(const Rcpp::NumericVector& x,
                                      const Rcpp::NumericVector& y,
                                      const Rcpp::NumericMatrix& z,
                                      const Rcpp::List& nb,
                                      const Rcpp::IntegerVector& libsizes,
                                      const Rcpp::IntegerVector& lib,
                                      const Rcpp::IntegerVector& pred,
                                      const Rcpp::IntegerVector& E,
                                      const Rcpp::IntegerVector& tau,
                                      const Rcpp::IntegerVector& b,
                                      bool simplex,
                                      double theta,
                                      int threads,
                                      int parallel_level,
                                      bool cumulate,
                                      bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> z_std(z.ncol());
  for (int i = 0; i < z.ncol(); ++i) {
    Rcpp::NumericVector covvar = z.column(i);
    z_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);
  std::vector<int> lib_std = Rcpp::as<std::vector<int>>(lib);
  std::vector<int> pred_std = Rcpp::as<std::vector<int>>(pred);
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Perform SCPCM For Lattice
  std::vector<std::vector<double>> result = SCPCM4Lattice(
    x_std,
    y_std,
    z_std,
    nb_vec,
    libsizes_std,
    lib_std,
    pred_std,
    E_std,
    tau_std,
    b_std,
    simplex,
    theta,
    threads,
    parallel_level,
    cumulate,
    progressbar);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 9);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
    resultMatrix(i, 5) = result[i][5];
    resultMatrix(i, 6) = result[i][6];
    resultMatrix(i, 7) = result[i][7];
    resultMatrix(i, 8) = result[i][8];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create(
    "libsizes","T_mean","D_mean",
    "T_sig","T_upper","T_lower",
    "D_sig","D_upper","D_lower");
  return resultMatrix;
}

// Wrapper function to perform GCMC Lattice and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCMC4Lattice(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::List& nb,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    const Rcpp::IntegerVector& b,
    const Rcpp::IntegerVector& max_r,
    int threads,
    bool progressbar){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_std = Rcpp::as<std::vector<int>>(lib);
  std::vector<int> pred_std = Rcpp::as<std::vector<int>>(pred);
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);
  std::vector<int> maxr_std = Rcpp::as<std::vector<int>>(max_r);

  // Remove values in b_std that are greater than validSampleNum or less than or equal to 3
  int validSampleNum = x_std.size();
  b_std.erase(std::remove_if(b_std.begin(), b_std.end(),
                             [validSampleNum](int x) { return x > validSampleNum || x <= 3; }),
                             b_std.end());

  // Remove duplicates for b_std
  std::sort(b_std.begin(), b_std.end());
  b_std.erase(std::unique(b_std.begin(), b_std.end()), b_std.end());

  // convert R based 1 index to C++ based 0 index
  for (size_t i = 0; i < lib_std.size(); ++i) {
    lib_std[i] -= 1;
  }
  for (size_t i = 0; i < pred_std.size(); ++i) {
    pred_std[i] -= 1;
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenLatticeEmbeddings(x_std, nb_vec, E[0], tau_std[0]);
  std::vector<std::vector<double>> e2 = GenLatticeEmbeddings(y_std, nb_vec, E[1], tau_std[1]);

  // Perform GCMC For Lattice
  std::vector<std::vector<double>> cs1 = CrossMappingCardinality(e1,e2,lib_std,pred_std,b_std,maxr_std,threads,progressbar);

  Rcpp::NumericMatrix resultMatrix(b_std.size(), 5);
  for (size_t i = 0; i < b_std.size(); ++i) {
    for (size_t j = 0; j < cs1[0].size(); ++j){
      resultMatrix(i, j) = cs1[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("neighbors",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_upper","x_xmap_y_lower");
  return resultMatrix;
}

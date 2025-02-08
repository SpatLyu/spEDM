#include <vector>
#include "CppLatticeUtils.h"
#include "Forecast4Lattice.h"
#include "GCCM4Lattice.h"
#include "SCPCM4Lattice.h"
#include "IntersectionCardinality.h"
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
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns: A NumericMatrix where each row contains {E, PearsonCor, MAE, RMSE}.
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSimplex4Lattice(const Rcpp::NumericVector& x,
                                        const Rcpp::List& nb,
                                        const Rcpp::IntegerVector& lib,
                                        const Rcpp::IntegerVector& pred,
                                        const Rcpp::IntegerVector& E,
                                        int b,
                                        int threads,
                                        bool includeself) {
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

  std::vector<std::vector<double>> res_std = Simplex4Lattice(
    vec_std,
    nb_vec,
    lib_indices,
    pred_indices,
    E_std,
    b,
    threads,
    includeself);

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

/**
 * Parameters:
 *   - x: A NumericVector containing the time series data.
 *   - nb: A List containing neighborhood information for lattice data.
 *   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
 *   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 *   - theta: A NumericVector containing the parameter values to be tested for theta.
 *   - E: An integer specifying the embedding dimension to test.
 *   - b: An integer specifying the number of neighbors to use for Simplex Mapping.
 *   - threads: An integer specifying the number of threads to use for parallel computation.
 *   - includeself: A boolean indicating whether to include the current state when constructing the embedding vector.
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
                                     int b,
                                     int threads,
                                     bool includeself) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

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

  std::vector<std::vector<double>> res_std = SMap4Lattice(
    vec_std,
    nb_vec,
    lib_indices,
    pred_indices,
    theta_std,
    E,
    b,
    threads,
    includeself);

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

// Wrapper function to perform GCCM Lattice and return a NumericMatrix
// predict y based on x ====> x xmap y ====> y causes x
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
                                     bool includeself,
                                     bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(x_std, nb_vec, E, includeself);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);

  // Define the interval [0, n-1] as a std::vector<std::pair<int, int>>
  int n = nb_vec.size();
  std::vector<std::pair<int, int>> interval = {{0, n-1}};

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
                                      const Rcpp::IntegerVector& E,
                                      int tau,
                                      int b,
                                      bool simplex,
                                      double theta,
                                      int threads,
                                      bool cumulate,
                                      bool includeself,
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
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  // Define the interval [0, n-1] as a std::vector<std::pair<int, int>>
  int n = nb_vec.size();
  std::vector<std::pair<int, int>> interval = {{0, n-1}};

  // Perform SCPCM For Lattice
  std::vector<std::vector<double>> result = SCPCM4Lattice(
    x_std,
    y_std,
    z_std,
    nb_vec,
    libsizes_std,
    interval,
    interval,
    E_std,
    tau,
    b,
    simplex,
    theta,
    threads,
    cumulate,
    includeself,
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
Rcpp::NumericVector RcppGCMC4Lattice(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::List& nb,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& b,
    const Rcpp::IntegerVector& max_r,
    int threads,
    bool includeself,
    bool progressbar){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);
  std::vector<int> maxr_std = Rcpp::as<std::vector<int>>(max_r);

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenLatticeEmbeddings(x_std, nb_vec, E[0], includeself);
  std::vector<std::vector<double>> e2 = GenLatticeEmbeddings(y_std, nb_vec, E[1], includeself);

  // Perform GCMC For Lattice
  double cs1 = IntersectionCardinality(e1,e2,b_std[0],maxr_std[0],threads,progressbar);
  double cs2 = IntersectionCardinality(e2,e1,b_std[1],maxr_std[1],threads,progressbar);

  Rcpp::NumericVector res_vec = Rcpp::NumericVector::create(
    Rcpp::Named("x_xmap_y",cs1),
    Rcpp::Named("y_xmap_x",cs2));

  return res_vec;
}

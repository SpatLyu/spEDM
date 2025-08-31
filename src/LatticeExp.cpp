#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "Forecast4Lattice.h"
#include "MultiViewEmbedding.h"
#include "GCCM4Lattice.h"
#include "SCPCM4Lattice.h"
#include "CrossMappingCardinality.h"
#include "FalseNearestNeighbors.h"
#include "SLM4Lattice.h"
#include "SGC4Lattice.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to convert Rcpp::List to std::vector<std::vector<int>> (the `nb` object)
std::vector<std::vector<int>> nb2vec(const Rcpp::List& nb) {
  // Get the number of elements in the nb object
  int n = nb.size();

  // Create a std::vector<std::vector<int>> to store the result
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

// Wrapper function to calculate accumulated lagged neighbor indices for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppLaggedNeighbor4Lattice(const Rcpp::List& nb, int lagNum = 1) {
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

// Wrapper function to calculate lagged values for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppLaggedVal4Lattice(const Rcpp::NumericVector& vec,
                                 const Rcpp::List& nb, int lagNum = 1) {
  int n = nb.size();

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Calculate lagged indices
  std::vector<std::vector<double>> lagged_values = CppLaggedVal4Lattice(vec_std, nb_vec, lagNum);

  // Convert std::vector<std::vector<int>> to Rcpp::List
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = Rcpp::wrap(lagged_values[i]);
  }

  return result;
}

// Wrapper function to generate embeddings for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenLatticeEmbeddings(const Rcpp::NumericVector& vec,
                                             const Rcpp::List& nb,
                                             int E = 3,
                                             int tau = 1,
                                             int style = 1) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec_std, nb_vec, E, tau, style);

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

// Wrapper function to generate neighbors for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGenLatticeNeighbors(const Rcpp::NumericVector& vec,
                                   const Rcpp::List& nb,
                                   const Rcpp::IntegerVector& lib,
                                   int k = 8) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector lib to std::vector<int> lib_std
  std::vector<int> lib_std;

  // Check that lib indices are within bounds & convert R based 1 index to C++ based 0 index
  int nsample = vec_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > nsample) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }

    if (!std::isnan(vec_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }

  // Generate neighbors
  std::vector<std::vector<int>> neighbors = GenLatticeNeighbors(
    vec_std, nb_vec, lib_std, static_cast<size_t>(std::abs(k))
  );

  // Convert neighbors to Rcpp::List with 1-based indexing
  int n = neighbors.size();
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int> neighbor_i = neighbors[i];
    for (auto& idx : neighbor_i) {
      idx += 1; // convert to 1-based index
    }
    result[i] = Rcpp::IntegerVector(neighbor_i.begin(), neighbor_i.end());
  }

  return result;
}

// Wrapper function to implement a symbolic transformation of a univariate spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenLatticeSymbolization(const Rcpp::NumericVector& vec,
                                                const Rcpp::List& nb,
                                                const Rcpp::IntegerVector& lib,
                                                const Rcpp::IntegerVector& pred,
                                                int k = 8) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_std;
  std::vector<int> pred_std;

  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  int n = vec_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(vec_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(vec_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  //  Generate symbolization map
  std::vector<double> symbolmap = GenLatticeSymbolization(
    vec_std, nb_vec, lib_std, pred_std, static_cast<size_t>(std::abs(k))
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(symbolmap);
}

// Wrapper function to partition spatial units in spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppDivideLattice(const Rcpp::List& nb,int b = 3) {
  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  //  Divide a spatial lattice into connected blocks
  std::vector<int> blocks = CppDivideLattice(nb_vec, b);

  // Convert the result back to Rcpp::IntegerVector
  return Rcpp::wrap(blocks);
}

// Wrapper function to perform univariate Spatial Logistic Map for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSLMUni4Lattice(
    const Rcpp::NumericVector& vec,
    const Rcpp::List& nb,
    int k = 4,
    int step = 20,
    double alpha = 0.77,
    double escape_threshold = 1e10
) {
  // Convert vec to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Call the core function
  std::vector<std::vector<double>> result = SLMUni4Lattice(vec_std, nb_vec, k, step, alpha, escape_threshold);

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result.size());
  int n_cols = step + 1;
  Rcpp::NumericMatrix out(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out(i, j) = result[i][j];
    }
  }

  return out;
}

// Wrapper function to perform bivariate Spatial Logistic Map for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMBi4Lattice(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::List& nb,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double beta_xy = 0.05,
    double beta_yx = 0.4,
    int interact = 0,
    double escape_threshold = 1e10
) {
  // Convert x/y to std::vector<double>
  std::vector<double> vec1 = Rcpp::as<std::vector<double>>(x);
  std::vector<double> vec2 = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMBi4Lattice(
    vec1, vec2, nb_vec, k, step, alpha_x, alpha_y, beta_xy, beta_yx, interact, escape_threshold
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = step + 1;
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y
  );

  return out;
}

// Wrapper function to perform trivariate Spatial Logistic Map for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMTri4Lattice(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& z,
    const Rcpp::List& nb,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double alpha_z = 0.55,
    double beta_xy = 0.05,
    double beta_xz = 0.05,
    double beta_yx = 0.4,
    double beta_yz = 0.4,
    double beta_zx = 0.65,
    double beta_zy = 0.65,
    int interact = 0,
    double escape_threshold = 1e10
) {
  // Convert x/y to std::vector<double>
  std::vector<double> vec1 = Rcpp::as<std::vector<double>>(x);
  std::vector<double> vec2 = Rcpp::as<std::vector<double>>(y);
  std::vector<double> vec3 = Rcpp::as<std::vector<double>>(z);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMTri4Lattice(
    vec1, vec2, vec3, nb_vec,
    k, step, alpha_x, alpha_y, alpha_z,
    beta_xy, beta_xz, beta_yx, beta_yz, beta_zx, beta_zy,
    interact, escape_threshold
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = step + 1;
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);
  Rcpp::NumericMatrix out_z(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
      out_z(i, j) = result[2][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y,
    Rcpp::Named("z") = out_z
  );

  return out;
}

// Wrapper function to perform FNN for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppFNN4Lattice(
    const Rcpp::NumericVector& vec,
    const Rcpp::List& nb,
    const Rcpp::NumericVector& rt,
    const Rcpp::NumericVector& eps,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    int tau = 1,
    int style = 1,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp *Vector to std::vector<*>
  std::vector<double> rt_std = Rcpp::as<std::vector<double>>(rt);
  std::vector<double> eps_std = Rcpp::as<std::vector<double>>(eps);
  std::vector<int> lib_std;
  std::vector<int> pred_std;

  int validSampleNum = vec_std.size();
  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > validSampleNum) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(vec_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > validSampleNum) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(vec_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  // Generate embeddings
  std::vector<double> E_std = Rcpp::as<std::vector<double>>(E);
  int max_E = CppMax(E_std, true);
  std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec_std, nb_vec, max_E, tau, style);

  // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
  bool L1norm = (dist_metric == 1);

  // Perform FNN for spatial lattice data
  std::vector<double> fnn = CppFNN(embeddings,lib_std,pred_std,rt_std,eps_std,L1norm,threads,parallel_level);

  // Convert the result back to Rcpp::NumericVector and set names as "E:1", "E:2", ..., "E:n"
  Rcpp::NumericVector result = Rcpp::wrap(fnn);
  Rcpp::CharacterVector resnames(result.size());
  for (int i = 0; i < result.size(); ++i) {
    resnames[i] = "E:" + std::to_string(i + 1);
  }
  result.names() = resnames;

  return result;
}

/**
 * Description:
 *   Computes Simplex projection for lattice data and returns a matrix containing
 *   the embedding dimension (E), Pearson correlation coefficient (PearsonCor),
 *   mean absolute error (MAE), and root mean squared error (RMSE).
 *
 * Parameters:
 *   - source: A NumericVector containing the source spatial cross-sectional data to be embedded.
 *   - target: A NumericVector containing the source spatial cross-sectional data to be predicted.
 *   - nb: A List containing neighborhood information for lattice data.
 *   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
 *   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 *   - E: An IntegerVector specifying the embedding dimensions to test.
 *   - b: An IntegerVector specifying the numbers of neighbors to use for simplex projection.
 *   - tau: An integer specifying the step of spatial lags for prediction. Default is 1.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it).  Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
 *
 * Returns:
 *   A NumericMatrix where each row contains {E, b, PearsonCor, MAE, RMSE}:
 *   - E: The tested embedding dimension.
 *   - b: The tested numbers of neighbors
 *   - PearsonCor: The Pearson correlation coefficient between the predicted and actual values.
 *   - MAE: The mean absolute error between the predicted and actual values.
 *   - RMSE: The root mean squared error between the predicted and actual values.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSimplex4Lattice(const Rcpp::NumericVector& source,
                                        const Rcpp::NumericVector& target,
                                        const Rcpp::List& nb,
                                        const Rcpp::IntegerVector& lib,
                                        const Rcpp::IntegerVector& pred,
                                        const Rcpp::IntegerVector& E,
                                        const Rcpp::IntegerVector& b,
                                        int tau = 1,
                                        int style = 1,
                                        int dist_metric = 2,
                                        bool dist_average = true,
                                        int threads = 8) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) && !std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) && !std::isnan(target_std[pred[i] - 1])) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
  }

  std::vector<std::vector<double>> res_std = Simplex4Lattice(
    source_std,
    target_std,
    nb_vec,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    style,
    dist_metric,
    dist_average,
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
 *   - source: A NumericVector containing the source spatial cross-sectional data to be embedded.
 *   - target: A NumericVector containing the source spatial cross-sectional data to be predicted.
 *   - nb: A List containing neighborhood information for lattice data.
 *   - lib: An IntegerVector specifying the library indices (1-based in R, converted to 0-based in C++).
 *   - pred: An IntegerVector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 *   - theta: A NumericVector containing the parameter values to be tested for theta.
 *   - E: The embedding dimension to evaluate. Default is 3.
 *   - tau: The spatial lag step for constructing lagged state-space vectors. Default is 1.
 *   - b: Number of nearest neighbors to use for prediction. Default is 4.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it). Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
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
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSMap4Lattice(const Rcpp::NumericVector& source,
                                     const Rcpp::NumericVector& target,
                                     const Rcpp::List& nb,
                                     const Rcpp::IntegerVector& lib,
                                     const Rcpp::IntegerVector& pred,
                                     const Rcpp::NumericVector& theta,
                                     int E = 3,
                                     int tau = 1,
                                     int b = 5,
                                     int style = 1,
                                     int dist_metric = 2,
                                     bool dist_average = true,
                                     int threads = 8) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) && !std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) && !std::isnan(target_std[pred[i] - 1])) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
  }

  std::vector<std::vector<double>> res_std = SMap4Lattice(
    source_std,
    target_std,
    nb_vec,
    lib_indices,
    pred_indices,
    theta_std,
    E,
    tau,
    b,
    style,
    dist_metric,
    dist_average,
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
 * - style: Embedding style selector (0: includes current state, 1: excludes it).
 * - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * - dist_average: Whether to average distance by the number of valid vector components.
 * - threads: An integer indicating the number of threads for parallel processing.
 *
 * Returns:
 * - An Rcpp::NumericVector containing the prediction results based on the multiview embedding.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppMultiView4Lattice(const Rcpp::NumericMatrix& x,
                                          const Rcpp::NumericVector& y,
                                          const Rcpp::List& nb,
                                          const Rcpp::IntegerVector& lib,
                                          const Rcpp::IntegerVector& pred,
                                          int E = 3,
                                          int tau = 1,
                                          int b = 5,
                                          int top = 5,
                                          int nvar = 3,
                                          int style = 1,
                                          int dist_metric = 2,
                                          int dist_average = true,
                                          int threads = 8){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> target = Rcpp::as<std::vector<double>>(y);

  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Initialize lib_indices and pred_indices with all false
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(target[lib[i] - 1])) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(target[pred[i] - 1])) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
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
    std::vector<std::vector<double>> vectors = GenLatticeEmbeddings(univec,nb_vec,E,tau,style);

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
    dist_metric,
    dist_average,
    threads);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(res);
}

// Wrapper function to compute intersection cardinality for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppIC4Lattice(const Rcpp::NumericVector& source,
                                   const Rcpp::NumericVector& target,
                                   const Rcpp::List& nb,
                                   const Rcpp::IntegerVector& lib,
                                   const Rcpp::IntegerVector& pred,
                                   const Rcpp::IntegerVector& E,
                                   const Rcpp::IntegerVector& b,
                                   int tau = 1,
                                   int exclude = 0,
                                   int style = 1,
                                   int dist_metric = 2,
                                   int threads = 8,
                                   int parallel_level = 0) {
  // Convert neighborhood list to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  std::vector<size_t> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) && !std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(static_cast<size_t>(lib[i] - 1)); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) && !std::isnan(target_std[pred[i] - 1])) {
      pred_indices.push_back(static_cast<size_t>(pred[i] - 1)); // Convert to 0-based index
    }
  }

  // Check the validity of the neignbor numbers
  std::vector<int> b_std;
  for (int i = 0; i < b.size(); ++i){
    if (b[i] > static_cast<int>(lib_indices.size())) {
      Rcpp::stop("Neighbor numbers count out of acceptable range at position %d (value: %d)", i + 1, b[i]);
    }
    b_std.push_back(b[i]);
  }

  std::vector<std::vector<double>> res_std = IC4Lattice(
    source_std,
    target_std,
    nb_vec,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    exclude,
    style,
    dist_metric,
    threads,
    parallel_level);

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
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "CausalScore", "Significance");
  return result;
}

// Wrapper function to perform GCCM for spatial lattice data
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGCCM4Lattice(const Rcpp::NumericVector& x,
                                     const Rcpp::NumericVector& y,
                                     const Rcpp::List& nb,
                                     const Rcpp::IntegerVector& libsizes,
                                     const Rcpp::IntegerVector& lib,
                                     const Rcpp::IntegerVector& pred,
                                     int E = 3,
                                     int tau = 1,
                                     int b = 5,
                                     bool simplex = true,
                                     double theta = 0,
                                     int threads = 8,
                                     int parallel_level = 0,
                                     int style = 1,
                                     int dist_metric = 2,
                                     bool dist_average = true,
                                     bool progressbar = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);
  std::vector<int> lib_std;
  std::vector<int> pred_std;

  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

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
    style,
    dist_metric,
    dist_average,
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

// Wrapper function to perform SCPCM for spatial lattice data
// predict y based on x ====> x xmap y ====> y causes x (account for controls)
// [[Rcpp::export(rng = false)]]
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
                                      bool simplex = true,
                                      double theta = 0,
                                      int threads = 8,
                                      int parallel_level = 0,
                                      bool cumulate = false,
                                      int style = 1,
                                      int dist_metric = 2,
                                      bool dist_average = true,
                                      bool single_sig = true,
                                      bool progressbar = false) {
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
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  std::vector<int> lib_std;
  std::vector<int> pred_std;
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

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
    style,
    dist_metric,
    dist_average,
    single_sig,
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

// Wrapper function to perform GCMC for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGCMC4Lattice(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::List& nb,
    const Rcpp::IntegerVector& libsizes,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int b = 4,
    int r = 0,
    int style = 1,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<size_t> libsizes_std = Rcpp::as<std::vector<size_t>>(libsizes);
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  int validSampleNum = x_std.size();
  // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  std::vector<size_t> lib_std;
  std::vector<size_t> pred_std;
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > validSampleNum) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > validSampleNum) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }
  }

  // check b that are greater than validSampleNum or less than or equal to 3
  if (b < 3 || b > validSampleNum) {
    Rcpp::stop("k cannot be less than or equal to 3 or greater than the number of non-NA values.");
  } else if (b + 1 > static_cast<int>(lib_std.size())){
    Rcpp::stop("Please check `libsizes` or `lib`; no valid libraries available for running GCMC.");
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenLatticeEmbeddings(x_std, nb_vec, E[0], tau_std[0], style);
  std::vector<std::vector<double>> e2 = GenLatticeEmbeddings(y_std, nb_vec, E[1], tau_std[1], style);

  // Perform GCMC for spatial lattice data
  CMCRes res = CrossMappingCardinality(e1,e2,libsizes_std,lib_std,pred_std,
                                       static_cast<size_t>(b),static_cast<size_t>(r),
                                       dist_metric,threads,parallel_level,progressbar);

  // Convert mean_aucs to Rcpp::DataFrame
  std::vector<double> libs, aucs;
  for (const auto& cm : res.causal_strength) {
    libs.push_back(cm[0]);
    aucs.push_back(cm[1]);
  }

  Rcpp::DataFrame cs_df = Rcpp::DataFrame::create(
    Rcpp::Named("libsizes") = libs,
    Rcpp::Named("x_xmap_y_mean") = aucs
  );

  // Wrap causal_strength with names
  Rcpp::DataFrame xmap_df = Rcpp::DataFrame::create(
    Rcpp::Named("neighbors") = res.cross_mapping[0],
    Rcpp::Named("x_xmap_y_mean") = res.cross_mapping[1],
    Rcpp::Named("x_xmap_y_sig") = res.cross_mapping[2],
    Rcpp::Named("x_xmap_y_upper") = res.cross_mapping[3],
    Rcpp::Named("x_xmap_y_lower")  = res.cross_mapping[4]
  );

  return Rcpp::List::create(
    Rcpp::Named("xmap") = xmap_df,
    Rcpp::Named("cs") = cs_df
  );
}

// Wrapper function to perform SGC for spatial lattice data without bootstrapped significance
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGCSingle4Lattice(const Rcpp::NumericVector& x,
                                          const Rcpp::NumericVector& y,
                                          const Rcpp::List& nb,
                                          const Rcpp::IntegerVector& lib,
                                          const Rcpp::IntegerVector& pred,
                                          int k,
                                          double base = 2,
                                          bool symbolize = true,
                                          bool normalize = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_std;
  std::vector<int> pred_std;

  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  // Perform SGC for spatial lattice data
  std::vector<double> sc = SGCSingle4Lattice(
    x_std,
    y_std,
    nb_vec,
    lib_std,
    pred_std,
    k,
    base,
    symbolize,
    normalize
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "statistic for y → x causality"
  );

  return sc_res;
}

// Wrapper function to perform SGC for spatial lattice data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGC4Lattice(const Rcpp::NumericVector& x,
                                    const Rcpp::NumericVector& y,
                                    const Rcpp::List& nb,
                                    const Rcpp::IntegerVector& lib,
                                    const Rcpp::IntegerVector& pred,
                                    const Rcpp::IntegerVector& block,
                                    int k,
                                    int threads = 8,
                                    int boot = 399,
                                    double base = 2,
                                    unsigned int seed = 42,
                                    bool symbolize = true,
                                    bool normalize = false,
                                    bool progressbar = true){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::List to std::vector<std::vector<int>>
  std::vector<std::vector<int>> nb_vec = nb2vec(nb);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_std;
  std::vector<int> pred_std;
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(block);

  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  // Perform SGC for spatial lattice data
  std::vector<double> sc = SGC4Lattice(
    x_std,
    y_std,
    nb_vec,
    lib_std,
    pred_std,
    b_std,
    k,
    threads,
    boot,
    base,
    seed,
    symbolize,
    normalize,
    progressbar
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "significance for x → y causality",
    "statistic for y → x causality",
    "significance for y → x causality"
   );

  return sc_res;
}

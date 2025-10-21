#include <vector>
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Computes predictions using the simplex projection method based on state-space reconstruction.
 *
 * See https://github.com/SpatLyu/simplex-smap-tutorial/blob/master/SimplexSmapFuncs.R for
 * the pure R Implementation
 *
 * Parameters:
 *   - embedding: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Spatial cross sectional series used as the target (should align with embedding).
 *   - lib: Integer vector of indices (which states to include when searching for neighbors, 1-based indexing).
 *   - pred: Integer vector of indices (which states to predict from, 1-based indexing).
 *   - num_neighbors: Number of neighbors to be used for simplex projection.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSimplexForecast(
    const Rcpp::NumericMatrix& embedding,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const int& dist_metric = 2,
    const bool& dist_average = true
  ){
  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> embedding_std(embedding.nrow(),
                                                 std::vector<double>(embedding.ncol()));
  for (int i = 0; i < embedding.nrow(); ++i) {
    for (int j = 0; j < embedding.ncol(); ++j) {
      embedding_std[i][j] = embedding(i, j);
    }
  }

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
  }

  // Call the SimplexProjectionPrediction function
  std::vector<double> pred_res = SimplexProjectionPrediction(
    embedding_std,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors,
    dist_metric,
    dist_average
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

/*
 * Computes predictions using the simplex projection method
 * for a collection of embeddings (multi-group combined version).
 *
 * This function extends RcppSimplexForecast by allowing the input
 * embeddings to be a list of matrices. Each matrix represents
 * one reconstructed state-space (rows = states, columns = embedding dimensions).
 * Internally, the list of matrices is converted to a 3-level nested
 * std::vector<std::vector<std::vector<double>>> and passed to
 * SimplexProjectionPrediction().
 *
 * Parameters:
 *   - embeddings: A list of R matrices, where each matrix represents
 *                 one reconstructed state-space embedding.
 *   - target: Numeric vector of target values (must align with embeddings).
 *   - lib: Integer vector (1-based indices) specifying library points.
 *   - pred: Integer vector (1-based indices) specifying prediction points.
 *   - num_neighbors: Number of nearest neighbors for simplex projection.
 *   - dist_metric: Distance metric (1 = Manhattan, 2 = Euclidean).
 *   - dist_average: Whether to average distance by valid vector components.
 *
 * Returns:
 *   Rcpp::NumericVector of predicted values (same length as target),
 *   with NaN entries for unpredicted indices.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSimplexForecastCom(
    const Rcpp::List& embeddings,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const int& dist_metric = 2,
    const bool& dist_average = true
) {
  // --- Convert embeddings (R list of matrices) into 3D std::vector
  std::vector<std::vector<std::vector<double>>> embeddings_std;
  embeddings_std.reserve(embeddings.size());

  for (int k = 0; k < embeddings.size(); ++k) {
    Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(embeddings[k]);
    std::vector<std::vector<double>> one_matrix(mat.nrow(), std::vector<double>(mat.ncol()));
    for (int i = 0; i < mat.nrow(); ++i) {
      for (int j = 0; j < mat.ncol(); ++j) {
        one_matrix[i][j] = mat(i, j);
      }
    }
    embeddings_std.push_back(std::move(one_matrix));
  }

  // --- Convert other arguments
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);
  std::vector<int> lib_indices, pred_indices;
  int n = target_std.size();

  // Convert R 1-based indices to C++ 0-based
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n)
      Rcpp::stop("lib index out of range: position %d (value %d)", i + 1, lib[i]);
    lib_indices.push_back(lib[i] - 1);
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n)
      Rcpp::stop("pred index out of range: position %d (value %d)", i + 1, pred[i]);
    pred_indices.push_back(pred[i] - 1);
  }

  // --- Call the C++ simplex projection function for composite embeddings
  std::vector<double> pred_res = SimplexProjectionPrediction(
    embeddings_std,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors,
    dist_metric,
    dist_average
  );

  return Rcpp::wrap(pred_res);
}

/*
 * Computes the S-Map forecast.
 *
 * See https://github.com/SpatLyu/simplex-smap-tutorial/blob/master/SimplexSmapFuncs.R for
 * the pure R Implementation
 *
 * Parameters:
 *   - embedding: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Spatial cross sectional series to be used as the target (should align with embedding).
 *   - lib: Integer vector of indices (which states to include when searching for neighbors, 1-based indexing).
 *   - pred: Integer vector of indices (which states to predict from, 1-based indexing).
 *   - num_neighbors: Number of neighbors to be used for S-Mapping.
 *   - theta: Weighting parameter for distances.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSMapForecast(
    const Rcpp::NumericMatrix& embedding,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const double& theta = 1.0,
    const int& dist_metric = 2,
    const bool& dist_average = true){
  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> embedding_std(embedding.nrow(),
                                                 std::vector<double>(embedding.ncol()));
  for (int i = 0; i < embedding.nrow(); ++i) {
    for (int j = 0; j < embedding.ncol(); ++j) {
      embedding_std[i][j] = embedding(i, j);
    }
  }

  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
  }

  // Call the SMapPrediction function
  std::vector<double> pred_res = SMapPrediction(
    embedding_std,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors,
    theta,
    dist_metric,
    dist_average
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

/*
 * Computes predictions using the S-Map method
 * for a collection of embeddings (multi-group combined version).
 *
 * This function extends RcppSMapForecast by allowing the input
 * embeddings to be a list of matrices. Each matrix represents
 * one reconstructed state-space (rows = states, columns = embedding dimensions).
 * Internally, the list of matrices is converted to a 3-level nested
 * std::vector<std::vector<std::vector<double>>> and passed to
 * SMapPrediction().
 *
 * Parameters:
 *   - embeddings: A list of R matrices, where each matrix represents
 *                 one reconstructed state-space embedding.
 *   - target: Numeric vector of target values (must align with embeddings).
 *   - lib: Integer vector (1-based indices) specifying library points.
 *   - pred: Integer vector (1-based indices) specifying prediction points.
 *   - num_neighbors: Number of nearest neighbors for local linear regression.
 *   - theta: Weighting parameter controlling exponential decay of distances.
 *   - dist_metric: Distance metric (1 = Manhattan, 2 = Euclidean).
 *   - dist_average: Whether to average distance by valid vector components.
 *
 * Returns:
 *   Rcpp::NumericVector of predicted values (same length as target),
 *   with NaN entries for unpredicted indices.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSMapForecastCom(
    const Rcpp::List& embeddings,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const double& theta = 1.0,
    const int& dist_metric = 2,
    const bool& dist_average = true
) {
  // --- Convert embeddings (R list of matrices) into 3D std::vector
  std::vector<std::vector<std::vector<double>>> embeddings_std;
  embeddings_std.reserve(embeddings.size());

  for (int k = 0; k < embeddings.size(); ++k) {
    Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(embeddings[k]);
    std::vector<std::vector<double>> one_matrix(mat.nrow(), std::vector<double>(mat.ncol()));
    for (int i = 0; i < mat.nrow(); ++i) {
      for (int j = 0; j < mat.ncol(); ++j) {
        one_matrix[i][j] = mat(i, j);
      }
    }
    embeddings_std.push_back(std::move(one_matrix));
  }

  // --- Convert other arguments
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);
  std::vector<int> lib_indices, pred_indices;
  int n = target_std.size();

  // Convert R 1-based indices to C++ 0-based
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n)
      Rcpp::stop("lib index out of range: position %d (value %d)", i + 1, lib[i]);
    lib_indices.push_back(lib[i] - 1);
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n)
      Rcpp::stop("pred index out of range: position %d (value %d)", i + 1, pred[i]);
    pred_indices.push_back(pred[i] - 1);
  }

  // --- Call the C++ S-Map prediction function for composite embeddings
  std::vector<double> pred_res = SMapPrediction(
    embeddings_std,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors,
    theta,
    dist_metric,
    dist_average
  );

  return Rcpp::wrap(pred_res);
}

/*
 * Computes the Intersection Cardinality (IC) curve
 *
 * This function serves as an interface between R and C++ to compute the Intersection Cardinality (IC) curve,
 * which quantifies the causal relationship between two variables by comparing the intersection of their nearest
 * neighbors in a state-space reconstruction. The function works by performing cross-mapping and calculating the
 * ratio of shared neighbors for each prediction index.
 *
 * Parameters:
 *   embedding_x: A NumericMatrix representing the state-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: A NumericMatrix representing the state-space reconstruction (embedded) of the potential effect variable.
 *   lib: An IntegerVector containing the library indices. These are 1-based indices in R, and will be converted to 0-based indices in C++.
 *   pred: An IntegerVector containing the prediction indices. These are 1-based indices in R, and will be converted to 0-based indices in C++.
 *   num_neighbors: An integer specifying the number of neighbors to use for cross mapping.
 *   n_excluded: An integer indicating the number of neighbors to exclude from the distance matrix.
 *   dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   threads: The number of parallel threads to use for computation.
 *   parallel_level: Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A NumericVector containing the intersection cardinality curve.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppIntersectionCardinality(
    const Rcpp::NumericMatrix& embedding_x,
    const Rcpp::NumericMatrix& embedding_y,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const int& n_excluded = 0,
    const int& dist_metric = 2,
    const int& threads = 8,
    const int& parallel_level = 0){
  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> e1(embedding_x.nrow(),
                                      std::vector<double>(embedding_x.ncol()));
  for (int i = 0; i < embedding_x.nrow(); ++i) {
    for (int j = 0; j < embedding_x.ncol(); ++j) {
      e1[i][j] = embedding_x(i, j);
    }
  }
  std::vector<std::vector<double>> e2(embedding_y.nrow(),
                                      std::vector<double>(embedding_y.ncol()));
  for (int i = 0; i < embedding_y.nrow(); ++i) {
    for (int j = 0; j < embedding_y.ncol(); ++j) {
      e2[i][j] = embedding_y(i, j);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  std::vector<size_t> pred_indices;

  int target_len = embedding_x.nrow();
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    lib_indices.push_back(static_cast<size_t>(lib[i] - 1)); // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    pred_indices.push_back(static_cast<size_t>(pred[i] - 1)); // Convert to 0-based index
  }

  if (lib_indices.size() < static_cast<size_t>(num_neighbors)){
    Rcpp::stop("Library size must not exceed the number of nearest neighbors used for mapping.");
  }

  // Call the IntersectionCardinality function
  std::vector<double> pred_res = IntersectionCardinality(
    e1,
    e2,
    lib_indices,
    pred_indices,
    static_cast<size_t>(num_neighbors),
    static_cast<size_t>(n_excluded),
    dist_metric,
    threads,
    parallel_level
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

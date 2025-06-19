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
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSimplexForecast(
    const Rcpp::NumericMatrix& embedding,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors){
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
    num_neighbors
  );

  // Convert the result back to Rcpp::NumericVector
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
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSMapForecast(
    const Rcpp::NumericMatrix& embedding,
    const Rcpp::NumericVector& target,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors,
    const double& theta){
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
    theta
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

/*
 * Computes the Intersection Cardinality (IC) scores
 *
 * This function serves as an interface between R and C++ to compute the Intersection Cardinality (IC) score,
 * which quantifies the causal relationship between two variables by comparing the intersection
 * of their nearest neighbors in a state-space reconstruction. The function works by performing cross-mapping
 * and calculating the ratio of shared neighbors for each prediction index.
 *
 * Parameters:
 *   embedding_x: A NumericMatrix representing the state-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: A NumericMatrix representing the state-space reconstruction (embedded) of the potential effect variable.
 *   lib: An IntegerVector containing the library indices. These are 1-based indices in R, and will be converted to 0-based indices in C++.
 *   pred: An IntegerVector containing the prediction indices. These are 1-based indices in R, and will be converted to 0-based indices in C++.
 *   num_neighbors: An integer specifying the number of neighbors to use for cross mapping.
 *   n_excluded: An integer indicating the number of neighbors to exclude from the distance matrix.
 *   threads: The number of parallel threads to use for computation.
 *   progressbar: A boolean value specifying whether to display a progress bar during computation.
 *
 * Returns:
 *   A NumericVector containing the intersection cardinality scores.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppIntersectionCardinality(
    const Rcpp::NumericMatrix& embedding_x,
    const Rcpp::NumericMatrix& embedding_y,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors,
    const int& n_excluded,
    const int& threads,
    const bool& progressbar){
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
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = embedding_x.nrow();
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

  // Call the IntersectionCardinality function
  std::vector<double> pred_res = IntersectionCardinality(
    e1,
    e2,
    lib_indices,
    pred_indices,
    num_neighbors,
    n_excluded,
    threads,
    progressbar
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

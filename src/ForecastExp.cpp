#include <vector>
#include "SimplexProjection.h"
#include "SMap.h"
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
 *   - target: Spatial cross-section series used as the target (should align with embedding).
 *   - lib: Integer vector of indices (which states to include when searching for neighbors, 1-based indexing).
 *   - pred: Integer vector of indices (which states to predict from, 1-based indexing).
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export]]
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

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(target_std.size(), false);
  std::vector<bool> pred_indices(target_std.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int libsize_int = lib.size();
  for (int i = 0; i < libsize_int; ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  int predsize_int = pred.size();
  for (int i = 0; i < predsize_int; ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
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
 *   - target: Spatial cross-section series to be used as the target (should align with embedding).
 *   - lib: Integer vector of indices (which states to include when searching for neighbors, 1-based indexing).
 *   - pred: Integer vector of indices (which states to predict from, 1-based indexing).
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: A Rcpp::NumericVector containing the predicted target values.
 */
// [[Rcpp::export]]
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

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(target_std.size(), false);
  std::vector<bool> pred_indices(target_std.size(), false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int libsize_int = lib.size();
  for (int i = 0; i < libsize_int; ++i) {
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  int predsize_int = pred.size();
  for (int i = 0; i < predsize_int; ++i) {
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
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

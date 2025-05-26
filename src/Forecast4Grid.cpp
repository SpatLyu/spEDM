#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for grid data using simplex projection.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbors to use for prediction.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given embedding dimension.
 */
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& source,
                                              const std::vector<std::vector<double>>& target,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<int>& E,
                                              const std::vector<int>& b,
                                              int tau,
                                              int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  int numRows = target.size();
  int numCols = target[0].size();

  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols); // Reserve space for efficiency

  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      vec_std.push_back(target[i][j]); // Add element to the vector
    }
  }

  // Sort and remove duplicates
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate unique pairs of E and b
  std::vector<std::pair<int, int>> unique_Ebcom;
  for (size_t i = 0; i < Es.size(); ++i){
    for (size_t j = 0; j < bs.size(); ++j){
      unique_Ebcom.emplace_back(Es[i],bs[j]); // Pair with E and b
    }
  }

  // Initialize result matrix with unique_Ebcom.size() rows and 5 columns
  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(5));

  // Parallel loop over each embedding dimension E
  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t i) {
    // Generate embeddings for the current E
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, unique_Ebcom[i].first, tau);

    // Compute metrics using SimplexBehavior
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, unique_Ebcom[i].second);

    // Store results in the matrix (no mutex needed since each thread writes to a unique index)
    result[i][0] = unique_Ebcom[i].first;   // Embedding dimension
    result[i][1] = unique_Ebcom[i].second;  // Number of nearest neighbors
    result[i][2] = metrics[0];              // Pearson correlation (rho)
    result[i][3] = metrics[1];              // Mean Absolute Error (MAE)
    result[i][4] = metrics[2];              // Root Mean Squared Error (RMSE)
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for grid data using the S-mapping method.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& source,
                                           const std::vector<std::vector<double>>& target,
                                           const std::vector<bool>& lib_indices,
                                           const std::vector<bool>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E,
                                           int tau,
                                           int b,
                                           int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  int numRows = target.size();
  int numCols = target[0].size();

  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols); // Reserve space for efficiency

  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      vec_std.push_back(target[i][j]); // Add element to the vector
    }
  }

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, E, tau);

  // Initialize result matrix with theta.size() rows and 4 columns
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  // Parallel loop over each theta parameter
  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {

    // Compute metrics using SimplexBehavior
    std::vector<double> metrics = SMapBehavior(embeddings, vec_std, lib_indices, pred_indices, b, theta[i]);

    // Store results in the matrix (no mutex needed since each thread writes to a unique index)
    result[i][0] = theta[i];           // Weighting parameter for distances
    result[i][1] = metrics[0];         // Pearson correlation (rho)
    result[i][2] = metrics[1];         // Mean Absolute Error (MAE)
    result[i][3] = metrics[2];         // Root Mean Squared Error (RMSE)
  }, threads_sizet);

  return result;
}

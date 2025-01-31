#include <vector>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

// Calculate the optimal embedding dimension of grid data using simplex projection
// Parameters:
//   - mat: A matrix will be embedded
//   - lib_indices: A boolean vector indicating library (training) set indices
//   - pred_indices: A boolean vector indicating prediction set indices
//   - E: A vector of embedding dimensions to evaluate
//   - b: Number of nearest neighbors to use for prediction
//   - threads: Number of threads used from the global pool
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns:
//   - A 2D vector where each row contains [E, rho, mae, rmse] for a given embedding dimension
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& mat,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<int>& E,
                                              double b,
                                              int threads,
                                              bool includeself) {
  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  int numRows = mat.size();
  int numCols = mat[0].size();

  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols); // Reserve space for efficiency

  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      vec_std.push_back(mat[i][j]); // Add element to the vector
    }
  }

  // Initialize result matrix with E.size() rows and 4 columns
  std::vector<std::vector<double>> result(E.size(), std::vector<double>(4));

  // Parallel loop over each embedding dimension E
  RcppThread::parallelFor(0, E.size(), [&](size_t i) {
    // Generate embeddings for the current E
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(mat, E[i], includeself);

    // Compute metrics using SimplexBehavior
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, b);

    // Store results in the matrix (no mutex needed since each thread writes to a unique index)
    result[i][0] = E[i];               // Embedding dimension
    result[i][1] = metrics[0];         // Pearson correlation (rho)
    result[i][2] = metrics[1];         // Mean Absolute Error (MAE)
    result[i][3] = metrics[2];         // Root Mean Squared Error (RMSE)
  }, threads_sizet);

  return result;
}

// Calculate the optimal theta parameter of grid data using s-mapping method
// Parameters:
//   - mat: A matrix will be embedded
//   - lib_indices: A boolean vector indicating library (training) set indices
//   - pred_indices: A boolean vector indicating prediction set indices
//   - theta: A vector of weighting parameters for distance calculation in SMap
//   - E: The embedding dimension to evaluate
//   - b: Number of nearest neighbors to use for prediction
//   - threads: Number of threads used from the global pool
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns:
//   - A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& mat,
                                           const std::vector<bool>& lib_indices,
                                           const std::vector<bool>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E,
                                           double b,
                                           int threads,
                                           bool includeself) {
  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  int numRows = mat.size();
  int numCols = mat[0].size();

  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols); // Reserve space for efficiency

  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      vec_std.push_back(mat[i][j]); // Add element to the vector
    }
  }

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(mat, E, includeself);

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

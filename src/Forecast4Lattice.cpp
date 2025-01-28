#include <vector>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

// Calculate the optimal embedding dimension of lattice data using simplex projection
// Parameters:
//   - vec: A vector of embedding values
//   - nb_vec: A 2D vector of neighbor indices
//   - lib_indices: A boolean vector indicating library (training) set indices
//   - pred_indices: A boolean vector indicating prediction set indices
//   - E: A vector of embedding dimensions to evaluate
//   - b: A parameter for the simplex behavior calculation
//   - threads: Number of threads used from the global pool
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns:
//   - A 2D vector where each row contains [E, rho, mae, rmse] for a given embedding dimension
std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& vec,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<bool>& lib_indices,
                                                 const std::vector<bool>& pred_indices,
                                                 const std::vector<int>& E,
                                                 double b,
                                                 int threads,
                                                 bool includeself) {
  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  // Initialize result matrix with E.size() rows and 4 columns
  std::vector<std::vector<double>> result(E.size(), std::vector<double>(4));

  // Parallel loop over each embedding dimension E
  RcppThread::parallelFor(0, E.size(), [&](size_t i) {
    // Generate embeddings for the current E
    std::vector<std::vector<double>> embeddings = GenLatticeEmbeddings(vec, nb_vec, E[i], includeself);

    // Compute metrics using SimplexBehavior
    std::vector<double> metrics = SimplexBehavior(embeddings, vec, lib_indices, pred_indices, b);

    // Store results in the matrix (no mutex needed since each thread writes to a unique index)
    result[i][0] = E[i];               // Embedding dimension
    result[i][1] = metrics[0];         // Pearson correlation (rho)
    result[i][2] = metrics[1];         // Mean Absolute Error (MAE)
    result[i][3] = metrics[2];         // Root Mean Squared Error (RMSE)
  }, threads_sizet);

  return result;
}

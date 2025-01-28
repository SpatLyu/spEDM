#include <vector>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

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

#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

/**
 * Calculate the one-dimensional index for a specified position in a 2D grid.
 *
 * This function converts a row and column position in a 2D grid into a corresponding
 * one-dimensional index, assuming the grid is stored in row-major order (i.e., rows are stored sequentially).
 *
 * @param curRow   The current row number (1-based indexing).
 * @param curCol   The current column number (1-based indexing).
 * @param totalRow The total number of rows in the grid.
 * @param totalCol The total number of columns in the grid.
 * @return         The calculated one-dimensional index (0-based indexing).
 */
int LocateGridIndices(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
}

/**
 * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size.
 *
 * This function calculates the cross mapping between a predictor variable (xEmbedings) and a response variable (yPred)
 * over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings   A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
 * @param yPred        A 1D vector of the response variable's values (spatial cross-section data).
 * @param lib_size     The size of the library (number of spatial units) used for prediction.
 * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param totalRow     The total number of rows in the 2D grid.
 * @param totalCol     The total number of columns in the 2D grid.
 * @param b            The number of nearest neighbors to use for prediction.
 * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @return             A vector of pairs, where each pair contains the library size and the corresponding cross mapping result.
 */
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta) {

  std::vector<std::pair<int, double>> x_xmap_y;
  double rho;

  for (int r = 1; r <= totalRow - lib_size + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size + 1; ++c) {

      // Initialize prediction and library indices
      std::vector<bool> pred_indices(totalRow * totalCol, false);
      std::vector<bool> lib_indices(totalRow * totalCol, false);

      // Set prediction indices
      for (const auto& p : pred) {
        pred_indices[LocateGridIndices(p.first, p.second, totalRow, totalCol)] = true;
      }

      // Exclude NA values in yPred from prediction indices
      for (size_t i = 0; i < yPred.size(); ++i) {
        if (std::isnan(yPred[i])) {
          pred_indices[i] = false;
        }
      }

      // Set library indices
      for (int i = r; i < r + lib_size; ++i) {
        for (int j = c; j < c + lib_size; ++j) {
          lib_indices[LocateGridIndices(i, j, totalRow, totalCol)] = true;
        }
      }

      // Check if more than half of the library is NA
      int na_count = 0;
      for (size_t i = 0; i < lib_indices.size(); ++i) {
        if (lib_indices[i] && std::isnan(yPred[i])) {
          ++na_count;
        }
      }
      if (na_count > (lib_size * lib_size) / 2) {
        continue;
      }

      // Run cross map and store results
      if (simplex){
        rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b);
      } else {
        rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta);
      }
      x_xmap_y.emplace_back(lib_size, rho);
    }
  }

  return x_xmap_y;
}

/**
 * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for multiple library sizes.
 *
 * This function calculates the cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
 * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
 *
 * @param xMatrix      A 2D matrix of the predictor variable's values (spatial cross-section data).
 * @param yMatrix      A 2D matrix of the response variable's values (spatial cross-section data).
 * @param lib_sizes    A vector of library sizes (number of spatial units) to use for prediction.
 * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param E            The number of dimensions for attractor reconstruction.
 * @param tau          The step of spatial lags for prediction.
 * @param b            The number of nearest neighbors to use for prediction.
 * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads      The number of threads to use for parallel processing.
 * @param includeself  Whether to include the current state when constructing the embedding vector.
 * @param progressbar  If true, display a progress bar during computation.
 * @return             A 2D vector where each row contains the library size, mean cross mapping result,
 *                     significance, and confidence interval bounds.
 */
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix, // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix, // Two dimension matrix of Y variable
    const std::vector<int>& lib_sizes,               // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& pred,    // Indices of spatial units to be predicted
    int E,                                           // Number of dimensions for the attractor reconstruction
    int tau,                                         // Step of spatial lags
    int b,                                           // Number of nearest neighbors to use for prediction
    bool simplex,                                    // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                    // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                     // Number of threads used from the global pool
    bool includeself,                                // Whether to include the current state when constructing the embedding vector
    bool progressbar                                 // Whether to print the progress bar
) {
  // If b is not provided correctly, default it to E + 2
  if (b <= 0) {
    b = E + 2;
  }

  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, E, includeself);

  // Ensure the maximum value does not exceed totalRow or totalCol
  int max_lib_size = std::max(totalRow, totalCol);

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is E + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, E + 2); });

  // Remove duplicates
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // // Iterate over each library size
  // for (int lib_size : unique_lib_sizes) {
  //   // Perform single grid cross-mapping for the current library size
  //   auto results = GCCMSingle4Grid(xEmbedings, yPred, lib_size, pred, totalRow, totalCol, b);
  //
  //   // Append the results to the main result container
  //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
  // }

  // Perform the operations using RcppThread
  if (progressbar) {
    RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = GCCMSingle4Grid(xEmbedings, yPred, lib_size, pred, totalRow, totalCol, b, simplex, theta);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = GCCMSingle4Grid(xEmbedings, yPred, lib_size, pred, totalRow, totalCol, b, simplex, theta);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    }, threads_sizet);
  }

  // Group by the first int and compute the mean
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  std::vector<std::vector<double>> final_results;
  for (const auto& group : grouped_results) {
    double mean_value = CppMean(group.second, true);
    final_results.push_back({static_cast<double>(group.first), mean_value});
  }

  int n = pred.size();
  // Calculate significance and confidence interval for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho = final_results[i][1];
    double significance = CppSignificance(rho, n);
    std::vector<double> confidence_interval = CppConfidence(rho, n);

    final_results[i].push_back(significance);
    final_results[i].push_back(confidence_interval[0]);
    final_results[i].push_back(confidence_interval[1]);
  }

  return final_results;
}

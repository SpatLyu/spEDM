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
 * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size and pred indice.
 *
 * This function calculates the cross mapping between a predictor variable (xEmbedings) and a response variable (yPred)
 * over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings     A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
 * @param yPred          A 1D vector of the response variable's values (spatial cross-section data).
 * @param lib_sizes      A vector of two integers, where the first element is the row-wise library size and the second element is the column-wise library size.
 * @param pred_indices   A boolean vector indicating which spatial units to be predicted.
 * @param totalRow       The total number of rows in the 2D grid.
 * @param totalCol       The total number of columns in the 2D grid.
 * @param b              The number of nearest neighbors to use for prediction.
 * @param simplex        If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta          The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads        The number of threads to use for parallel processing.
 * @param row_size_mark  If true, use the row-wise libsize to mark the libsize; if false, use col-wise libsize.
 *
 * @return  A vector of pairs, where each pair contains the library size and the corresponding cross mapping result.
 */
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& pred_indices,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    bool row_size_mark) {
  // Extract row-wise and column-wise library sizes
  int lib_size_row = lib_sizes[0];
  int lib_size_col = lib_sizes[1];

  // Determine the marked libsize
  int libsize = (row_size_mark) ? lib_size_row : lib_size_col;

  // Precompute valid (r, c) pairs
  std::vector<std::pair<int, int>> valid_indices;
  for (int r = 1; r <= totalRow - lib_size_row + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size_col + 1; ++c) {
      valid_indices.emplace_back(r, c);
    }
  }

  // // Initialize the result container with the same size as valid_indices
  // std::vector<std::pair<int, double>> x_xmap_y;
  // x_xmap_y.resize(valid_indices.size());

  // Preallocate the result vector to avoid out-of-bounds access
  std::vector<std::pair<int, double>> x_xmap_y(valid_indices.size());

  double rho;

  // // Iterate through precomputed (r, c) pairs
  // for (size_t i = 0; i < valid_indices.size(); ++i) {
  //   int r = valid_indices[i].first;
  //   int c = valid_indices[i].second;
  //
  //   // Initialize library indices
  //   std::vector<bool> lib_indices(totalRow * totalCol, false);
  //
  //   // Set library indices
  //   for (int i = r; i < r + lib_size_row; ++i) {
  //     for (int j = c; j < c + lib_size_col; ++j) {
  //       lib_indices[LocateGridIndices(i, j, totalRow, totalCol)] = true;
  //     }
  //   }
  //
  //   // Check if more than half of the library is NA
  //   int na_count = 0;
  //   for (size_t i = 0; i < lib_indices.size(); ++i) {
  //     if (lib_indices[i] && std::isnan(yPred[i])) {
  //       ++na_count;
  //     }
  //   }
  //
  //   if (na_count > (lib_size_row * lib_size_col) / 2) {
  //     rho = std::numeric_limits<int>::min();
  //   } else {
  //     // Run cross map and store results
  //     if (simplex) {
  //       rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b);
  //     } else {
  //       rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta);
  //     }
  //   }
  //
  //   std::pair<int, double> result(libsize, rho); // Store the product of row and column library sizes
  //   x_xmap_y[i] = result;
  // }

  // Perform the operations using RcppThread
  RcppThread::parallelFor(0, valid_indices.size(), [&](size_t i) {
    int r = valid_indices[i].first;
    int c = valid_indices[i].second;

    // Initialize library indices
    std::vector<bool> lib_indices(totalRow * totalCol, false);

    // Set library indices
    for (int i = r; i < r + lib_size_row; ++i) {
      for (int j = c; j < c + lib_size_col; ++j) {
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

    if (na_count > (lib_size_row * lib_size_col) / 2) {
      rho = std::numeric_limits<int>::min();
    } else {
      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b);
      } else {
        rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta);
      }
    }

    std::pair<int, double> result(libsize, rho); // Store the product of row and column library sizes
    x_xmap_y[i] = result;
  }, threads);

  return x_xmap_y;
}

/**
 * Perform Geographical Convergent Cross Mapping (GCCM) for spatial grid data.
 *
 * This function calculates the cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
 * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
 *
 * @param xMatrix      A 2D matrix of the predictor variable's values (spatial cross-section data).
 * @param yMatrix      A 2D matrix of the response variable's values (spatial cross-section data).
 * @param lib_sizes    A 2D vector where the first sub-vector contains row-wise library sizes and the second sub-vector contains column-wise library sizes.
 * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param E            The number of dimensions for attractor reconstruction.
 * @param tau          The step of spatial lags for prediction.
 * @param b            The number of nearest neighbors to use for prediction.
 * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads      The number of threads to use for parallel processing.
 * @param progressbar  If true, display a progress bar during computation.
 *
 * @return A 2D vector where each row contains the library size, mean cross mapping result,
 *         significance, and confidence interval bounds.
 */
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<std::vector<int>>& lib_sizes,
    const std::vector<std::pair<int, int>>& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    bool progressbar
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
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, E, tau);

  // Ensure the maximum value does not exceed totalRow or totalCol
  int max_lib_size_row = totalRow;
  int max_lib_size_col = totalCol;

  // Extract row-wise and column-wise library sizes
  std::vector<int> row_lib_sizes = lib_sizes[0];
  std::vector<int> col_lib_sizes = lib_sizes[1];

  // Transform to ensure no size exceeds max_lib_size_row or max_lib_size_col
  std::transform(row_lib_sizes.begin(), row_lib_sizes.end(), row_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_row); });
  std::transform(col_lib_sizes.begin(), col_lib_sizes.end(), col_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_col); });

  // Remove duplicates in row-wise and column-wise library sizes
  row_lib_sizes.erase(std::unique(row_lib_sizes.begin(), row_lib_sizes.end()), row_lib_sizes.end());
  col_lib_sizes.erase(std::unique(col_lib_sizes.begin(), col_lib_sizes.end()), col_lib_sizes.end());

  // // Generate unique pairs of row-wise and column-wise library sizes
  // std::vector<std::pair<int, int>> unique_lib_size_pairs;
  // for (int row_size : row_lib_sizes) {
  //   for (int col_size : col_lib_sizes) {
  //     unique_lib_size_pairs.emplace_back(row_size, col_size);
  //   }
  // }

  // Generate unique pairs of row-wise and column-wise library sizes
  std::vector<std::pair<int, int>> unique_lib_size_pairs;

  // Determine which library size vector is longer
  int row_size_count = row_lib_sizes.size();
  int col_size_count = col_lib_sizes.size();
  int min_size = std::min(row_size_count, col_size_count);
  // int max_size = std::max(row_size_count, col_size_count);

  // Fill unique_lib_size_pairs based on the shorter vector
  for (int i = 0; i < min_size; ++i) {
    unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes[i]);
  }

  bool row_size_mark = true;
  // Handle the excess elements for the longer vector
  if (row_size_count > col_size_count) {
    for (int i = min_size; i < row_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes.back()); // Pair with the max value of col_lib_sizes
    }
  }

  if (row_size_count < col_size_count) {
    for (int i = min_size; i < col_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes.back(), col_lib_sizes[i]); // Pair with the max value of row_lib_sizes
    }
    row_size_mark = false;
  }

  // Set prediction indices
  std::vector<bool> pred_indices(totalRow * totalCol, false);
  for (const auto& p : pred) {
    pred_indices[LocateGridIndices(p.first, p.second, totalRow, totalCol)] = true;
  }

  // Exclude NA values in yPred from prediction indices
  for (size_t i = 0; i < yPred.size(); ++i) {
    if (std::isnan(yPred[i])) {
      pred_indices[i] = false;
    }
  }

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // Iterate over each library size
  if (progressbar) {
    RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
    for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
      int lib_size_row = unique_lib_size_pairs[i].first;
      int lib_size_col = unique_lib_size_pairs[i].second;
      auto results = GCCMSingle4Grid(xEmbedings, yPred, {lib_size_row, lib_size_col}, pred_indices, totalRow, totalCol, b, simplex, theta, threads_sizet, row_size_mark);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
      bar++;
    }
  } else {
    for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
      int lib_size_row = unique_lib_size_pairs[i].first;
      int lib_size_col = unique_lib_size_pairs[i].second;
      auto results = GCCMSingle4Grid(xEmbedings, yPred, {lib_size_row, lib_size_col}, pred_indices, totalRow, totalCol, b, simplex, theta, threads_sizet, row_size_mark);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    }
  }

  // Group by the first int (library size) and compute the mean
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
    double significance = CppCorSignificance(rho, n);
    std::vector<double> confidence_interval = CppCorConfidence(rho, n);

    final_results[i].push_back(significance);
    final_results[i].push_back(confidence_interval[0]);
    final_results[i].push_back(confidence_interval[1]);
  }

  return final_results;
}

// #include <vector>
// #include <algorithm>
// #include <cmath>
// #include <numeric>
// #include <utility>
// #include <limits>
// #include <map>
// #include "CppStats.h"
// #include "CppGridUtils.h"
// #include "SimplexProjection.h"
// #include "SMap.h"
// #include <RcppThread.h>
//
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(RcppThread)]]
//
// /**
//  * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size and pred indice.
//  *
//  * This function calculates the cross mapping between a predictor variable (xEmbedings) and a response variable (yPred)
//  * over a 2D grid, using either Simplex Projection or S-Mapping.
//  *
//  * @param xEmbedings     A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
//  * @param yPred          A 1D vector of the response variable's values (spatial cross-section data).
//  * @param lib_sizes      A vector of two integers, where the first element is the row-wise library size and the second element is the column-wise library size.
//  * @param pred_indices   A boolean vector indicating which spatial units to be predicted.
//  * @param totalRow       The total number of rows in the 2D grid.
//  * @param totalCol       The total number of columns in the 2D grid.
//  * @param b              The number of nearest neighbors to use for prediction.
//  * @param simplex        If true, use Simplex Projection; if false, use S-Mapping.
//  * @param theta          The distance weighting parameter for S-Mapping (ignored if simplex is true).
//  * @param row_size_mark  If ture, use the row-wise libsize to mark the libsize; if false, use col-wise libsize.
//  *
//  * @return  A vector of pairs, where each pair contains the library size and the corresponding cross mapping result.
//  */
// std::vector<std::pair<int, double>> GCCMSingle4Grid(
//     const std::vector<std::vector<double>>& xEmbedings,
//     const std::vector<double>& yPred,
//     const std::vector<int>& lib_sizes,
//     const std::vector<bool>& pred_indices,
//     int totalRow,
//     int totalCol,
//     int b,
//     bool simplex,
//     double theta,
//     bool row_size_mark) {
//
//   std::vector<std::pair<int, double>> x_xmap_y;
//   double rho;
//
//   // Extract row-wise and column-wise library sizes
//   int lib_size_row = lib_sizes[0];
//   int lib_size_col = lib_sizes[1];
//
//   // Determine the marked libsize
//   int libsize = (row_size_mark) ? lib_size_row : lib_size_col;
//
//   for (int r = 1; r <= totalRow - lib_size_row + 1; ++r) {
//     for (int c = 1; c <= totalCol - lib_size_col + 1; ++c) {
//
//       // Initialize library indices
//       std::vector<bool> lib_indices(totalRow * totalCol, false);
//
//       // Set library indices
//       for (int i = r; i < r + lib_size_row; ++i) {
//         for (int j = c; j < c + lib_size_col; ++j) {
//           lib_indices[LocateGridIndices(i, j, totalRow, totalCol)] = true;
//         }
//       }
//
//       // Check if more than half of the library is NA
//       int na_count = 0;
//       for (size_t i = 0; i < lib_indices.size(); ++i) {
//         if (lib_indices[i] && std::isnan(yPred[i])) {
//           ++na_count;
//         }
//       }
//       if (na_count > (lib_size_row * lib_size_col) / 2) {
//         continue;
//       }
//
//       // Run cross map and store results
//       if (simplex){
//         rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b);
//       } else {
//         rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta);
//       }
//       x_xmap_y.emplace_back(libsize, rho); // Store the product of row and column library sizes
//     }
//   }
//
//   return x_xmap_y;
// }
//
// /**
//  * Perform Geographical Convergent Cross Mapping (GCCM) for spatial grid data.
//  *
//  * This function calculates the cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
//  * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
//  *
//  * @param xMatrix      A 2D matrix of the predictor variable's values (spatial cross-section data).
//  * @param yMatrix      A 2D matrix of the response variable's values (spatial cross-section data).
//  * @param lib_sizes    A 2D vector where the first sub-vector contains row-wise library sizes and the second sub-vector contains column-wise library sizes.
//  * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
//  * @param E            The number of dimensions for attractor reconstruction.
//  * @param tau          The step of spatial lags for prediction.
//  * @param b            The number of nearest neighbors to use for prediction.
//  * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
//  * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
//  * @param threads      The number of threads to use for parallel processing.
//  * @param progressbar  If true, display a progress bar during computation.
//  *
//  * @return A 2D vector where each row contains the library size, mean cross mapping result,
//  *         significance, and confidence interval bounds.
//  */
// std::vector<std::vector<double>> GCCM4Grid(
//     const std::vector<std::vector<double>>& xMatrix,
//     const std::vector<std::vector<double>>& yMatrix,
//     const std::vector<std::vector<int>>& lib_sizes,
//     const std::vector<std::pair<int, int>>& pred,
//     int E,
//     int tau,
//     int b,
//     bool simplex,
//     double theta,
//     int threads,
//     bool progressbar
// ) {
//   // If b is not provided correctly, default it to E + 2
//   if (b <= 0) {
//     b = E + 2;
//   }
//
//   size_t threads_sizet = static_cast<size_t>(threads);
//   unsigned int max_threads = std::thread::hardware_concurrency();
//   threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);
//
//   // Get the dimensions of the xMatrix
//   int totalRow = xMatrix.size();
//   int totalCol = xMatrix[0].size();
//
//   // Flatten yMatrix into a 1D array (row-major order)
//   std::vector<double> yPred;
//   for (const auto& row : yMatrix) {
//     yPred.insert(yPred.end(), row.begin(), row.end());
//   }
//
//   // Generate embeddings for xMatrix
//   std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, E, tau);
//
//   // Ensure the maximum value does not exceed totalRow or totalCol
//   int max_lib_size_row = totalRow;
//   int max_lib_size_col = totalCol;
//
//   // Extract row-wise and column-wise library sizes
//   std::vector<int> row_lib_sizes = lib_sizes[0];
//   std::vector<int> col_lib_sizes = lib_sizes[1];
//
//   // Transform to ensure no size exceeds max_lib_size_row or max_lib_size_col
//   std::transform(row_lib_sizes.begin(), row_lib_sizes.end(), row_lib_sizes.begin(),
//                  [&](int size) { return std::min(size, max_lib_size_row); });
//   std::transform(col_lib_sizes.begin(), col_lib_sizes.end(), col_lib_sizes.begin(),
//                  [&](int size) { return std::min(size, max_lib_size_col); });
//
//   // Remove duplicates in row-wise and column-wise library sizes
//   row_lib_sizes.erase(std::unique(row_lib_sizes.begin(), row_lib_sizes.end()), row_lib_sizes.end());
//   col_lib_sizes.erase(std::unique(col_lib_sizes.begin(), col_lib_sizes.end()), col_lib_sizes.end());
//
//   // // Generate unique pairs of row-wise and column-wise library sizes
//   // std::vector<std::pair<int, int>> unique_lib_size_pairs;
//   // for (int row_size : row_lib_sizes) {
//   //   for (int col_size : col_lib_sizes) {
//   //     unique_lib_size_pairs.emplace_back(row_size, col_size);
//   //   }
//   // }
//
//   // Generate unique pairs of row-wise and column-wise library sizes
//   std::vector<std::pair<int, int>> unique_lib_size_pairs;
//
//   // Determine which library size vector is longer
//   int row_size_count = row_lib_sizes.size();
//   int col_size_count = col_lib_sizes.size();
//   int min_size = std::min(row_size_count, col_size_count);
//   // int max_size = std::max(row_size_count, col_size_count);
//
//   // Fill unique_lib_size_pairs based on the shorter vector
//   for (int i = 0; i < min_size; ++i) {
//     unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes[i]);
//   }
//
//   bool row_size_mark = true;
//   // Handle the excess elements for the longer vector
//   if (row_size_count > col_size_count) {
//     for (int i = min_size; i < row_size_count; ++i) {
//       unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes.back()); // Pair with the max value of col_lib_sizes
//     }
//   }
//
//   if (row_size_count < col_size_count) {
//     for (int i = min_size; i < col_size_count; ++i) {
//       unique_lib_size_pairs.emplace_back(row_lib_sizes.back(), col_lib_sizes[i]); // Pair with the max value of row_lib_sizes
//     }
//     row_size_mark = false;
//   }
//
//   // Set prediction indices
//   std::vector<bool> pred_indices(totalRow * totalCol, false);
//   for (const auto& p : pred) {
//     pred_indices[LocateGridIndices(p.first, p.second, totalRow, totalCol)] = true;
//   }
//
//   // Exclude NA values in yPred from prediction indices
//   for (size_t i = 0; i < yPred.size(); ++i) {
//     if (std::isnan(yPred[i])) {
//       pred_indices[i] = false;
//     }
//   }
//
//   // Initialize the result container
//   std::vector<std::pair<int, double>> x_xmap_y;
//
//   // // Iterate over each library size
//   // for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
//   //   int lib_size_row = unique_lib_size_pairs[i].first;
//   //   int lib_size_col = unique_lib_size_pairs[i].second;
//   //   auto results = GCCMSingle4Grid(xEmbedings, yPred, {lib_size_row, lib_size_col}, pred_indices, totalRow, totalCol, b, simplex, theta, row_size_mark);
//   //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
//   // }
//
//   // Perform the operations using RcppThread
//   if (progressbar) {
//     RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
//     RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
//       int lib_size_row = unique_lib_size_pairs[i].first;
//       int lib_size_col = unique_lib_size_pairs[i].second;
//       auto results = GCCMSingle4Grid(xEmbedings, yPred, {lib_size_row, lib_size_col}, pred_indices, totalRow, totalCol, b, simplex, theta, row_size_mark);
//       x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
//       bar++;
//     }, threads_sizet);
//   } else {
//     RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
//       int lib_size_row = unique_lib_size_pairs[i].first;
//       int lib_size_col = unique_lib_size_pairs[i].second;
//       auto results = GCCMSingle4Grid(xEmbedings, yPred, {lib_size_row, lib_size_col}, pred_indices, totalRow, totalCol, b, simplex, theta, row_size_mark);
//       x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
//     }, threads_sizet);
//   }
//
//   // Group by the first int (library size) and compute the mean
//   std::map<int, std::vector<double>> grouped_results;
//   for (const auto& result : x_xmap_y) {
//     grouped_results[result.first].push_back(result.second);
//   }
//
//   std::vector<std::vector<double>> final_results;
//   for (const auto& group : grouped_results) {
//     double mean_value = CppMean(group.second, true);
//     final_results.push_back({static_cast<double>(group.first), mean_value});
//   }
//
//   int n = pred.size();
//   // Calculate significance and confidence interval for each result
//   for (size_t i = 0; i < final_results.size(); ++i) {
//     double rho = final_results[i][1];
//     double significance = CppCorSignificance(rho, n);
//     std::vector<double> confidence_interval = CppCorConfidence(rho, n);
//
//     final_results[i].push_back(significance);
//     final_results[i].push_back(confidence_interval[0]);
//     final_results[i].push_back(confidence_interval[1]);
//   }
//
//   return final_results;
// }

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
#include <RcppThread.h>

// Function to locate the index in a 2D grid
int locate(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
}

// GCCMSingle4Grid function
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b) {

  std::vector<std::pair<int, double>> x_xmap_y;

  for (int r = 1; r <= totalRow - lib_size + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size + 1; ++c) {

      // Initialize prediction and library indices
      std::vector<bool> pred_indices(totalRow * totalCol, false);
      std::vector<bool> lib_indices(totalRow * totalCol, false);

      // Set prediction indices
      for (const auto& p : pred) {
        pred_indices[locate(p.first, p.second, totalRow, totalCol)] = true;
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
          lib_indices[locate(i, j, totalRow, totalCol)] = true;
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
      double results = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b);
      x_xmap_y.emplace_back(lib_size, results);
    }
  }

  return x_xmap_y;
}

// GCCM4Grid function
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<int>& lib_sizes,
    const std::vector<std::pair<int, int>>& pred,
    int E,
    int tau = 1,
    int b = 0) {

  // If b is not provided, default it to E + 1
  if (b == 0) {
    b = E + 1;
  }

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, E);

  // Ensure the minimum value in unique_lib_sizes is E + 2 and remove duplicates
  std::vector<int> unique_lib_sizes;
  std::unique_copy(lib_sizes.begin(), lib_sizes.end(), std::back_inserter(unique_lib_sizes),
                   [&](int a, int b) { return a == b; });
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::max(size, E + 2); });

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
  RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
  RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
    int lib_size = unique_lib_sizes[i];
    auto results = GCCMSingle4Grid(xEmbedings, yPred, lib_size, pred, totalRow, totalCol, b);
    x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    bar++;
  });

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

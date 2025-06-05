#include <cmath>
#include <limits>
#include <vector>
#include <numeric>
#include "CppGridUtils.h"

/**
 * @brief Simulate a univariate Spatial Logistic Map (SLM) over grid-structured data.
 *
 * This function performs time-stepped simulations of the Spatial Logistic Map
 * over a 2D grid. Each cell evolves based on its own value and the average
 * of its k nearest Queen-adjacent neighbors (excluding NaNs).
 *
 * @param mat                2D grid of initial values (e.g., population densities), row-major.
 * @param k                  Number of neighbors to consider (using Queen adjacency).
 * @param step               Number of simulation time steps to run.
 * @param alpha              Growth/interaction parameter in the logistic update rule.
 * @param escape_threshold   Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A 2D vector of simulation results:
 *         Each row corresponds to a spatial unit (flattened from the grid),
 *         and each column to a time step (0 to step).
 */
std::vector<std::vector<double>> SLMUni4Grid(
    const std::vector<std::vector<double>>& mat,
    size_t k,
    size_t step,
    double alpha,
    double escape_threshold = 1e10
){
  size_t nrow = mat.size();
  size_t ncol = mat[0].size();
  size_t ncell = nrow * ncol;

  // Convert 2D grid to 1D vector
  std::vector<double> vec(ncell, std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < nrow; ++i){
    for (size_t j = 0; j < ncol; ++j){
      vec[i * ncol + j] = mat[i][j];
    }
  }

  // Initialize index library for all spatial units
  std::vector<int> lib(ncell);
  std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., ncell-1

  // Build k-nearest neighbors for each cell based on Queen adjacency
  std::vector<std::vector<int>> neighbors = GenGridNeighbors(mat, lib, k);

  // Initialize result matrix with NaNs
  std::vector<std::vector<double>> res(ncell,
                                       std::vector<double>(step + 1,
                                                           std::numeric_limits<double>::quiet_NaN()));

  // Set initial values
  for (size_t i = 0; i < ncell; ++i){
    res[i][0] = vec[i];
  }

  // Time-stepped simulation
  for (size_t s = 1; s <= step; ++s){
    for (size_t i = 0; i < ncell; ++i){
      if (std::isnan(res[i][s - 1])) continue;

      double v_neighbors = 0.0;
      double valid_neighbors = 0.0;
      const std::vector<int>& local_neighbors = neighbors[i];

      for (size_t j = 0; j < local_neighbors.size(); ++j){
        int neighbor_idx = local_neighbors[j];
        if (!std::isnan(res[neighbor_idx][s - 1])){
          v_neighbors += res[neighbor_idx][s - 1];
          valid_neighbors += 1.0;
        }
      }

      double v_next = std::numeric_limits<double>::quiet_NaN();
      if (valid_neighbors > 0){
        v_next = 1 - alpha * res[i][s - 1] * v_neighbors / valid_neighbors;
      }

      if (std::abs(v_next) <= escape_threshold){
        res[i][s] = v_next;
      }
    }
  }

  return res;
}

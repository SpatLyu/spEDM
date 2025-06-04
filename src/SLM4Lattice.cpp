#include <cmath>
#include <limits>
#include <vector>
#include <numeric>
#include "CppLatticeUtils.h"

/**
 * @brief Simulate a univariate Spatial Logistic Map (SLM) over lattice-structured data.
 *
 * This function performs time-stepped simulations of the Spatial Logistic Map
 * on a lattice data where each spatial unit evolves based on its own value and
 * the average of its k nearest neighbors.
 *
 * @param vec                Initial values of the lattice data (e.g., population densities).
 * @param nb                 Neighbor list for each lattice unit (e.g., rook or queen adjacency).
 * @param k                  Number of neighbors to consider.
 * @param step               Number of simulation time steps to run.
 * @param alpha              Growth/interaction parameter in the logistic update rule.
 * @param escape_threshold   Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A 2D vector of simulation results:
 *         Each row corresponds to a spatial unit,
 *         and each column to a time step (0 to step).
 */
std::vector<std::vector<double>> SLMUni4Lattice(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    size_t step,
    double alpha,
    double escape_threshold = 1e10
){
  // Initialize index library for all spatial units
  std::vector<int> lib(vec.size());
  // for(size_t i = 0; i < vec.size(); ++i){
  //   lib[i] = static_cast<int>(i);
  // }
  std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., vec.size()-1

  // Generate fixed-k neighbors (if possible)
  std::vector<std::vector<int>> neighbors = GenLatticeNeighbors(vec, nb, lib, k);

  // Initialize result matrix with NaNs (rows: spatial units, cols: time steps)
  std::vector<std::vector<double>> res(vec.size(),
                                       std::vector<double>(step + 1,
                                                           std::numeric_limits<double>::quiet_NaN()));

  // Set initial values at time step 0
  for(size_t i = 0; i < vec.size(); ++i){
    res[i][0] = vec[i];
  }

  // Time-stepped simulation
  for (size_t s = 1; s <= step; ++s){
    for(size_t currentIndex = 0; currentIndex < vec.size(); ++currentIndex){
      // Skip if the current value is invalid (NaN)
      if (std::isnan(res[currentIndex][s - 1])) continue;

      // Compute the average of valid neighboring values
      double v_neighbors = 0;
      double valid_neighbors = 0;
      const std::vector<int>& local_neighbors = neighbors[currentIndex];
      for (size_t i = 0; i < local_neighbors.size(); ++i) {
        if (!std::isnan(res[local_neighbors[i]][s - 1])){
          v_neighbors += res[local_neighbors[i]][s - 1];
          valid_neighbors += 1;
        }
      }

      // Apply the spatial logistic map update if neighbors exist
      double v_next = std::numeric_limits<double>::quiet_NaN();
      if (valid_neighbors > 0){
        v_next = 1 - alpha * res[currentIndex][s - 1] * v_neighbors / valid_neighbors;
      }

      // Update result only if the value is within the escape threshold
      if (std::abs(v_next) <= escape_threshold){
        res[currentIndex][s] = v_next;
      }
    }
  }

  return res;
}

/**
 * @brief Simulate a bivariate Spatial Logistic Map (SLM) over lattice-structured data.
 *
 * This function performs time-stepped simulations of a coupled bivariate Spatial Logistic Map
 * on lattice data, where each of the two spatial variables evolves based on its own previous value,
 * the average of its neighbors' values, and cross-variable interaction from the other variable.
 *
 * For each spatial unit, the evolution of variable 1 is influenced by its own spatial neighbors
 * and an inhibitory term proportional to variable 2 at the same location, and vice versa.
 *
 * @param vec1               Initial values of the first spatial variable (e.g., species A density).
 * @param vec2               Initial values of the second spatial variable (e.g., species B density).
 * @param nb                 Neighbor list for each spatial unit (e.g., rook or queen adjacency).
 * @param k                  Number of neighbors to consider.
 * @param step               Number of simulation time steps to run.
 * @param alpha1             Growth/interaction parameter for the first variable.
 * @param alpha2             Growth/interaction parameter for the second variable.
 * @param beta1              Cross-inhibition coefficient from variable 2 to variable 1.
 * @param beta2              Cross-inhibition coefficient from variable 1 to variable 2.
 * @param escape_threshold   Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A 3D vector of simulation results:
 *         - First dimension: variable index (0 for vec1, 1 for vec2),
 *         - Second dimension: spatial units,
 *         - Third dimension: time steps (0 to step).
 */
std::vector<std::vector<std::vector<double>>> SLMBi4Lattice(
    const std::vector<double>& vec1,
    const std::vector<double>& vec2,
    const std::vector<std::vector<int>>& nb,
    size_t k,
    size_t step,
    double alpha1,
    double alpha2,
    double beta1,
    double beta2,
    double escape_threshold = 1e10
){
  // Initialize index library for all spatial units
  std::vector<int> lib(vec1.size());
  // for(size_t i = 0; i < vec1.size(); ++i){
  //   lib[i] = static_cast<int>(i);
  // }
  std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., vec1.size()-1

  // Generate fixed-k neighbors (if possible)
  std::vector<std::vector<int>> neighbors = GenLatticeNeighbors(vec1, nb, lib, k);

  // Initialize result array with NaNs (2, rows: spatial units, cols: time steps)
  std::vector<std::vector<std::vector<double>>> res(2,
                                                    std::vector<std::vector<double>>(vec1.size(),
                                                                                     std::vector<double>(step + 1,
                                                                                                         std::numeric_limits<double>::quiet_NaN())));

  // Set initial values at time step 0
  for(size_t i = 0; i < vec1.size(); ++i){
    res[0][i][0] = vec1[i];
    res[1][i][0] = vec2[i];
  }

  // Time-stepped simulation
  for (size_t s = 1; s <= step; ++s){
    for(size_t currentIndex = 0; currentIndex < vec1.size(); ++currentIndex){
      // Skip if the current value is invalid (NaN)
      if (std::isnan(res[0][currentIndex][s - 1]) && std::isnan(res[1][currentIndex][s - 1])) continue;

      // Compute the average of valid neighboring values
      double v_neighbors_1 = 0;
      double v_neighbors_2 = 0;
      double valid_neighbors_1 = 0;
      double valid_neighbors_2 = 0;
      const std::vector<int>& local_neighbors = neighbors[currentIndex];
      for (size_t i = 0; i < local_neighbors.size(); ++i) {
        if (!std::isnan(res[0][local_neighbors[i]][s - 1])){
          v_neighbors_1 += res[0][local_neighbors[i]][s - 1];
          valid_neighbors_1 += 1;
        }
        if (!std::isnan(res[1][local_neighbors[i]][s - 1])){
          v_neighbors_2 += res[1][local_neighbors[i]][s - 1];
          valid_neighbors_2 += 1;
        }
      }

      // Apply the spatial logistic map update if neighbors exist
      double v_next_1 = std::numeric_limits<double>::quiet_NaN();
      double v_next_2 = std::numeric_limits<double>::quiet_NaN();
      if (valid_neighbors_1 > 0){
        v_next_1 = 1 - alpha1 * res[0][currentIndex][s - 1] * (v_neighbors_1 / valid_neighbors_1 - beta1 * res[1][currentIndex][s - 1]);
      }
      if (valid_neighbors_2 > 0){
        v_next_2 = 1 - alpha2 * res[1][currentIndex][s - 1] * (v_neighbors_2 / valid_neighbors_2 - beta2 * res[0][currentIndex][s - 1]);
      }

      // Update result only if the value is within the escape threshold
      if (std::abs(v_next_1) <= escape_threshold){
        res[0][currentIndex][s] = v_next_1;
      }
      if (std::abs(v_next_2) <= escape_threshold){
        res[1][currentIndex][s] = v_next_1;
      }
    }
  }

  return res;
}

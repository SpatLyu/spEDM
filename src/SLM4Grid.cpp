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

  // Initialize result matrix with NaNs
  std::vector<std::vector<double>> res(ncell,
                                       std::vector<double>(step + 1,
                                                           std::numeric_limits<double>::quiet_NaN()));

  // Set initial values
  for (size_t i = 0; i < ncell; ++i){
    res[i][0] = vec[i];
  }

  if (k > 0){
    // Initialize index library for all spatial units
    std::vector<int> lib(ncell);
    std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., ncell-1

    // Build k-nearest neighbors for each cell based on Queen adjacency
    std::vector<std::vector<int>> neighbors = GenGridNeighbors(mat, lib, k);

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
  } else {
    // Time-stepped simulation
    for (size_t s = 1; s <= step; ++s){
      for (size_t i = 0; i < ncell; ++i){
        if (std::isnan(res[i][s - 1])) continue;

        // Apply the logistic map update if no neighbors exist
        double v_next = res[i][s - 1] * (alpha - alpha * res[i][s - 1]);

        if (std::abs(v_next) <= escape_threshold){
          res[i][s] = v_next;
        }
      }
    }
  }

  return res;
}

/**
 * @brief Simulate a bivariate Spatial Logistic Map (SLM) on gridded data.
 *
 * This function performs time-stepped simulations of a bivariate Spatial Logistic Map
 * on two input grid-based datasets. Each spatial unit (cell) evolves over time based on
 * its own state and the spatial influence from its k nearest neighbors, determined via
 * Queen-style adjacency and expanded recursively until k non-NaN neighbors are obtained.
 *
 * The two interacting variables influence each other through cross-inhibition terms
 * (beta12 and beta21), enabling simulation of coupled spatial dynamics (e.g., predator-prey,
 * competing species, or interacting fields).
 *
 * The logistic update rule for each cell includes a local term and a spatial term,
 * and diverging values (e.g., explosions) are filtered using an escape threshold.
 *
 * @param mat1              Initial grid values for the first variable.
 * @param mat2              Initial grid values for the second variable.
 * @param k                 Number of valid spatial neighbors to use (Queen adjacency with recursive expansion).
 * @param step              Number of simulation time steps to run.
 * @param alpha1            Growth/interaction parameter for the first variable.
 * @param alpha2            Growth/interaction parameter for the second variable.
 * @param beta12            Cross-inhibition from the first variable to the second.
 * @param beta21            Cross-inhibition from the second variable to the first.
 * @param escape_threshold  Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A 3D vector of simulation results:
 *         - Dimension 0: variable index (0 for mat1, 1 for mat2),
 *         - Dimension 1: spatial unit index (row-major order),
 *         - Dimension 2: time step (0 to step).
 */
std::vector<std::vector<std::vector<double>>> SLMBi4Grid(
    const std::vector<std::vector<double>>& mat1,
    const std::vector<std::vector<double>>& mat2,
    size_t k,
    size_t step,
    double alpha1,
    double alpha2,
    double beta12,
    double beta21,
    double escape_threshold = 1e10
){
  size_t nrow = mat1.size();
  size_t ncol = mat1[0].size();
  size_t ncell = nrow * ncol;

  // Convert 2D grid to 1D vector
  std::vector<double> vec1(ncell, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> vec2(ncell, std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < nrow; ++i){
    for (size_t j = 0; j < ncol; ++j){
      vec1[i * ncol + j] = mat1[i][j];
      vec2[i * ncol + j] = mat2[i][j];
    }
  }

  // Initialize result array with NaNs (2, rows: spatial units, cols: time steps)
  std::vector<std::vector<std::vector<double>>> res(2,
                                                    std::vector<std::vector<double>>(ncell,
                                                                                     std::vector<double>(step + 1,
                                                                                                         std::numeric_limits<double>::quiet_NaN())));

  // Set initial values at time step 0
  for(size_t i = 0; i < vec1.size(); ++i){
    res[0][i][0] = vec1[i];
    res[1][i][0] = vec2[i];
  }

  if (k > 0){
    // Initialize index library for all spatial units
    std::vector<int> lib(ncell);
    std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., ncell-1

    // Build k-nearest neighbors for each cell based on Queen adjacency
    std::vector<std::vector<int>> neighbors = GenGridNeighbors(mat1, lib, k);

    // Time-stepped simulation
    for (size_t s = 1; s <= step; ++s){
      for (size_t i = 0; i < ncell; ++i){
        // Skip if the current value is invalid (NaN)
        if (std::isnan(res[0][i][s - 1]) && std::isnan(res[1][i][s - 1])) continue;

        // Compute the average of valid neighboring values
        double v_neighbors_1 = 0;
        double v_neighbors_2 = 0;
        double valid_neighbors_1 = 0;
        double valid_neighbors_2 = 0;
        const std::vector<int>& local_neighbors = neighbors[i];

        for (size_t j = 0; j < local_neighbors.size(); ++j){
          int neighbor_idx = local_neighbors[j];
          if (!std::isnan(res[0][neighbor_idx][s - 1])){
            v_neighbors_1 += res[0][neighbor_idx][s - 1];
            valid_neighbors_1 += 1.0;
          }
          if (!std::isnan(res[1][neighbor_idx][s - 1])){
            v_neighbors_2 += res[1][neighbor_idx][s - 1];
            valid_neighbors_2 += 1.0;
          }
        }

        // Apply the spatial logistic map update if neighbors exist
        double v_next_1 = std::numeric_limits<double>::quiet_NaN();
        double v_next_2 = std::numeric_limits<double>::quiet_NaN();
        if (valid_neighbors_1 > 0){
          v_next_1 = 1 - alpha1 * res[0][i][s - 1] * (v_neighbors_1 / valid_neighbors_1 - beta21 * res[1][i][s - 1]);
        }
        if (valid_neighbors_2 > 0){
          v_next_2 = 1 - alpha2 * res[1][i][s - 1] * (v_neighbors_2 / valid_neighbors_2 - beta12 * res[0][i][s - 1]);
        }

        // Update result only if the value is within the escape threshold
        if (std::abs(v_next_1) <= escape_threshold){
          res[0][i][s] = v_next_1;
        }
        if (std::abs(v_next_2) <= escape_threshold){
          res[1][i][s] = v_next_2;
        }
      }
    }
  } else {
    // Time-stepped simulation
    for (size_t s = 1; s <= step; ++s){
      for (size_t i = 0; i < ncell; ++i){
        // Skip if the current value is invalid (NaN)
        if (std::isnan(res[0][i][s - 1]) && std::isnan(res[1][i][s - 1])) continue;

        // Apply the logistic map update if no neighbors exist
        double v_next_1 = res[0][i][s - 1] * (alpha1 - alpha1 * res[0][i][s - 1] - beta21 * res[1][i][s - 1]);
        double v_next_2 = res[1][i][s - 1] * (alpha2 - alpha2 * res[1][i][s - 1] - beta12 * res[0][i][s - 1]);

        // Update result only if the value is within the escape threshold
        if (std::abs(v_next_1) <= escape_threshold){
          res[0][i][s] = v_next_1;
        }
        if (std::abs(v_next_2) <= escape_threshold){
          res[1][i][s] = v_next_2;
        }
      }
    }
  }

  return res;
}

/**
 * @brief Simulate a three-variable spatial logistic map on a 2D grid.
 *
 * This function performs a time-stepped simulation of three interacting spatial variables
 * arranged in grid format. For each cell in the spatial lattice, it constructs k-nearest neighbors
 * based on Queen adjacency. At each time step, the function updates each variable's value by applying
 * a coupled logistic map formula that considers the average values of neighboring cells and the influence
 * of the other two variables through specified interaction coefficients.
 *
 * The simulation proceeds for a specified number of steps, starting from initial input matrices,
 * and stops or skips updates when values become invalid (NaN) or exceed a defined escape threshold.
 *
 * @param mat1 Initial values for variable 1 in 2D grid form.
 * @param mat2 Initial values for variable 2 in 2D grid form.
 * @param mat3 Initial values for variable 3 in 2D grid form.
 * @param k Number of nearest neighbors to consider for spatial interaction.
 * @param step Number of time steps for the simulation.
 * @param alpha1 Growth rate parameter for variable 1.
 * @param alpha2 Growth rate parameter for variable 2.
 * @param alpha3 Growth rate parameter for variable 3.
 * @param beta12 Interaction coefficient from variable 1 to variable 2.
 * @param beta13 Interaction coefficient from variable 1 to variable 3.
 * @param beta21 Interaction coefficient from variable 2 to variable 1.
 * @param beta23 Interaction coefficient from variable 2 to variable 3.
 * @param beta31 Interaction coefficient from variable 3 to variable 1.
 * @param beta32 Interaction coefficient from variable 3 to variable 2.
 * @param escape_threshold Threshold to prevent values from diverging too far.
 *
 * @return A 3D vector containing simulated values for each variable,
 *         spatial unit, and time step.
 */
std::vector<std::vector<std::vector<double>>> SLMTri4Grid(
    const std::vector<std::vector<double>>& mat1,
    const std::vector<std::vector<double>>& mat2,
    const std::vector<std::vector<double>>& mat3,
    size_t k,
    size_t step,
    double alpha1,
    double alpha2,
    double alpha3,
    double beta12,
    double beta13,
    double beta21,
    double beta23,
    double beta31,
    double beta32,
    double escape_threshold = 1e10
){
  size_t nrow = mat1.size();
  size_t ncol = mat1[0].size();
  size_t ncell = nrow * ncol;

  // Convert 2D grid to 1D vector
  std::vector<double> vec1(ncell, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> vec2(ncell, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> vec3(ncell, std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < nrow; ++i){
    for (size_t j = 0; j < ncol; ++j){
      vec1[i * ncol + j] = mat1[i][j];
      vec2[i * ncol + j] = mat2[i][j];
      vec3[i * ncol + j] = mat3[i][j];
    }
  }

  // Initialize result array with NaNs (3, rows: spatial units, cols: time steps)
  std::vector<std::vector<std::vector<double>>> res(3,
                                                    std::vector<std::vector<double>>(ncell,
                                                                                     std::vector<double>(step + 1,
                                                                                                         std::numeric_limits<double>::quiet_NaN())));

  // Set initial values at time step 0
  for(size_t i = 0; i < vec1.size(); ++i){
    res[0][i][0] = vec1[i];
    res[1][i][0] = vec2[i];
    res[2][i][0] = vec3[i];
  }

  if (k > 0){
    // Initialize index library for all spatial units
    std::vector<int> lib(ncell);
    std::iota(lib.begin(), lib.end(), 0); // Fill with 0, 1, ..., ncell-1

    // Build k-nearest neighbors for each cell based on Queen adjacency
    std::vector<std::vector<int>> neighbors = GenGridNeighbors(mat1, lib, k);

    // Time-stepped simulation
    for (size_t s = 1; s <= step; ++s){
      for (size_t i = 0; i < ncell; ++i){
        // Skip if the current value is invalid (NaN)
        if (std::isnan(res[0][i][s - 1]) &&
            std::isnan(res[1][i][s - 1]) &&
            std::isnan(res[2][i][s - 1])) continue;

        // Compute the average of valid neighboring values
        double v_neighbors_1 = 0;
        double v_neighbors_2 = 0;
        double v_neighbors_3 = 0;
        double valid_neighbors_1 = 0;
        double valid_neighbors_2 = 0;
        double valid_neighbors_3 = 0;
        const std::vector<int>& local_neighbors = neighbors[i];

        for (size_t j = 0; j < local_neighbors.size(); ++j){
          int neighbor_idx = local_neighbors[j];
          if (!std::isnan(res[0][neighbor_idx][s - 1])){
            v_neighbors_1 += res[0][neighbor_idx][s - 1];
            valid_neighbors_1 += 1.0;
          }
          if (!std::isnan(res[1][neighbor_idx][s - 1])){
            v_neighbors_2 += res[1][neighbor_idx][s - 1];
            valid_neighbors_2 += 1.0;
          }
          if (!std::isnan(res[2][neighbor_idx][s - 1])){
            v_neighbors_3 += res[2][neighbor_idx][s - 1];
            valid_neighbors_3 += 1.0;
          }
        }

        // Apply the spatial logistic map update if neighbors exist
        double v_next_1 = std::numeric_limits<double>::quiet_NaN();
        double v_next_2 = std::numeric_limits<double>::quiet_NaN();
        double v_next_3 = std::numeric_limits<double>::quiet_NaN();
        if (valid_neighbors_1 > 0){
          v_next_1 = 1 - alpha1 * res[0][i][s - 1] * (v_neighbors_1 / valid_neighbors_1 - beta21 * res[1][i][s - 1] - beta31 * res[2][i][s - 1]);
        }
        if (valid_neighbors_2 > 0){
          v_next_2 = 1 - alpha2 * res[1][i][s - 1] * (v_neighbors_2 / valid_neighbors_2 - beta12 * res[0][i][s - 1] - beta32 * res[2][i][s - 1]);
        }
        if (valid_neighbors_3 > 0){
          v_next_3 = 1 - alpha3 * res[2][i][s - 1] * (v_neighbors_3 / valid_neighbors_3 - beta13 * res[0][i][s - 1] - beta23 * res[1][i][s - 1]);
        }

        // Update result only if the value is within the escape threshold
        if (std::abs(v_next_1) <= escape_threshold){
          res[0][i][s] = v_next_1;
        }
        if (std::abs(v_next_2) <= escape_threshold){
          res[1][i][s] = v_next_2;
        }
        if (std::abs(v_next_3) <= escape_threshold){
          res[2][i][s] = v_next_3;
        }
      }
    }
  } else {
    // Time-stepped simulation
    for (size_t s = 1; s <= step; ++s){
      for (size_t i = 0; i < ncell; ++i){
        // Skip if the current value is invalid (NaN)
        if (std::isnan(res[0][i][s - 1]) &&
            std::isnan(res[1][i][s - 1]) &&
            std::isnan(res[2][i][s - 1])) continue;

        // Apply the logistic map update if no neighbors exist
        double v_next_1 = res[0][i][s - 1] * (alpha1 - alpha1 * res[0][i][s - 1] - beta21 * res[1][i][s - 1] - beta31 * res[2][i][s - 1]);
        double v_next_2 = res[1][i][s - 1] * (alpha2 - alpha2 * res[1][i][s - 1] - beta12 * res[0][i][s - 1] - beta32 * res[2][i][s - 1]);
        double v_next_3 = res[2][i][s - 1] * (alpha3 - alpha3 * res[2][i][s - 1] - beta13 * res[0][i][s - 1] - beta23 * res[1][i][s - 1]);

        // Update result only if the value is within the escape threshold
        if (std::abs(v_next_1) <= escape_threshold){
          res[0][i][s] = v_next_1;
        }
        if (std::abs(v_next_2) <= escape_threshold){
          res[1][i][s] = v_next_2;
        }
        if (std::abs(v_next_3) <= escape_threshold){
          res[2][i][s] = v_next_3;
        }
      }
    }
  }

  return res;
}

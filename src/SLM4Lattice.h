#ifndef SLM4Lattice_H
#define SLM4Lattice_H

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
);

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
 * @param beta12             Cross-inhibition coefficient from variable 1 to variable 2.
 * @param beta21             Cross-inhibition coefficient from variable 2 to variable 1.
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
    double beta12,
    double beta21,
    double escape_threshold = 1e10
);

/**
 * @brief Simulate a trivariate Spatial Logistic Map (SLM) over lattice-structured data.
 *
 * This function simulates the dynamics of a three-variable coupled Spatial Logistic Map
 * across a lattice. Each spatial variable evolves over discrete time steps under the
 * influence of: (1) its own previous value, (2) the mean of its spatial neighbors,
 * and (3) cross-variable interactions from the other two variables at the same location.
 *
 * For each spatial unit:
 * - Variable 1 is influenced by its own neighbors and is inhibited by Variable 2 and Variable 3.
 * - Variable 2 is influenced by its own neighbors and is inhibited by Variable 1 and Variable 3.
 * - Variable 3 is influenced by its own neighbors and is inhibited by Variable 1 and Variable 2.
 *
 * @param vec1               Initial values of the first spatial variable.
 * @param vec2               Initial values of the second spatial variable.
 * @param vec3               Initial values of the third spatial variable.
 * @param nb                 Neighbor list for each spatial unit (e.g., rook or queen adjacency).
 * @param k                  Number of neighbors to select per unit (fixed-k).
 * @param step               Number of simulation time steps to perform.
 * @param alpha1             Growth/interaction parameter for variable 1.
 * @param alpha2             Growth/interaction parameter for variable 2.
 * @param alpha3             Growth/interaction parameter for variable 3.
 * @param beta12             Cross-inhibition from variable 1 to variable 2.
 * @param beta13             Cross-inhibition from variable 1 to variable 3.
 * @param beta21             Cross-inhibition from variable 2 to variable 1.
 * @param beta23             Cross-inhibition from variable 2 to variable 3.
 * @param beta31             Cross-inhibition from variable 3 to variable 1.
 * @param beta32             Cross-inhibition from variable 3 to variable 2.
 * @param escape_threshold   Threshold beyond which values are treated as divergent (default: 1e10).
 *
 * @return A 3D vector of simulation results:
 *         - First dimension: variable index (0 for vec1, 1 for vec2, 2 for vec3),
 *         - Second dimension: spatial units,
 *         - Third dimension: time steps (0 to step).
 */
std::vector<std::vector<std::vector<double>>> SLMTri4Lattice(
    const std::vector<double>& vec1,
    const std::vector<double>& vec2,
    const std::vector<double>& vec3,
    const std::vector<std::vector<int>>& nb,
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
);

#endif // SLM4Lattice_H

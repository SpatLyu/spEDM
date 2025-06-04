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

#endif // SLM4Lattice_H

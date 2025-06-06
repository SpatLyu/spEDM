#ifndef SLM4Grid_H
#define SLM4Grid_H

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
);

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
);

/**
 * @brief Simulate a three-variable spatial logistic map on a 2D grid.
 *
 * This function performs a time-stepped simulation of three interacting spatial variables
 * arranged in grid format. For each cell in the spatial grid, it constructs k-nearest neighbors
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
);

#endif // SLM4Grid_H

#ifndef Forecast4Lattice_H
#define Forecast4Lattice_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for lattice data using simplex projection.
 *
 * Parameters:
 *   - vec: A vector to be embedded.
 *   - nb_vec: A 2D vector of neighbor indices.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbor values to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given combination of embedding dimension and nearest neighbors.
 */
std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& vec,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<bool>& lib_indices,
                                                 const std::vector<bool>& pred_indices,
                                                 const std::vector<int>& E,
                                                 const std::vector<int>& b,
                                                 int tau,
                                                 int threads);

/*
 * Evaluates prediction performance of different theta parameters for lattice data using the s-mapping method.
 *
 * Parameters:
 *   - vec: A vector to be embedded.
 *   - nb_vec: A 2D vector of neighbor indices.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Lattice(const std::vector<double>& vec,
                                              const std::vector<std::vector<int>>& nb_vec,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<double>& theta,
                                              int E,
                                              int tau,
                                              int b,
                                              int threads);

#endif // Forecast4Lattice_H

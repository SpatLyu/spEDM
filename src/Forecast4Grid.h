#ifndef Forecast4Grid_H
#define Forecast4Grid_H

#include <vector>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different embedding dimensions for grid data using simplex projection.
 *
 * Parameters:
 *   - mat: A matrix to be embedded.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *   - includeself: Whether to include the current state when constructing the embedding vector.
 *
 * Returns:
 *   A 2D vector where each row contains [E, rho, mae, rmse] for a given embedding dimension.
 */
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& mat,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<int>& E,
                                              double b,
                                              int threads,
                                              bool includeself);

/*
 * Evaluates prediction performance of different theta parameters for grid data using the S-mapping method.
 *
 * Parameters:
 *   - mat: A matrix to be embedded.
 *   - lib_indices: A boolean vector indicating library (training) set indices.
 *   - pred_indices: A boolean vector indicating prediction set indices.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *   - includeself: Whether to include the current state when constructing the embedding vector.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& mat,
                                           const std::vector<bool>& lib_indices,
                                           const std::vector<bool>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E,
                                           double b,
                                           int threads,
                                           bool includeself);

#endif // Forecast4Grid_H

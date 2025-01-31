#ifndef Forecast4Lattice_H
#define Forecast4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// Calculate the optimal embedding dimension of lattice data using simplex projection
// Parameters:
//   - vec: A vector of embedding values
//   - nb_vec: A 2D vector of neighbor indices
//   - lib_indices: A boolean vector indicating library (training) set indices
//   - pred_indices: A boolean vector indicating prediction set indices
//   - E: A vector of embedding dimensions to evaluate
//   - b: Number of nearest neighbors to use for prediction
//   - threads: Number of threads used from the global pool
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns:
//   - A 2D vector where each row contains [E, rho, mae, rmse] for a given embedding dimension
std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& vec,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<bool>& lib_indices,
                                                 const std::vector<bool>& pred_indices,
                                                 const std::vector<int>& E,
                                                 double b,
                                                 int threads,
                                                 bool includeself);

// Calculate the optimal theta parameter of lattice data using s-mapping method
// Parameters:
//   - vec: A vector will be embedded
//   - nb_vec: A 2D vector of neighbor indices
//   - lib_indices: A boolean vector indicating library (training) set indices
//   - pred_indices: A boolean vector indicating prediction set indices
//   - theta: A vector of weighting parameters for distance calculation in SMap
//   - E: The embedding dimension to evaluate
//   - b: Number of nearest neighbors to use for prediction
//   - threads: Number of threads used from the global pool
//   - includeself: Whether to include the current state when constructing the embedding vector
// Returns:
//   - A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value
std::vector<std::vector<double>> SMap4Lattice(const std::vector<double>& vec,
                                              const std::vector<std::vector<int>>& nb_vec,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<double>& theta,
                                              int E,
                                              double b,
                                              int threads,
                                              bool includeself);

#endif // Forecast4Lattice_H

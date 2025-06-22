#ifndef Forecast4Grid_H
#define Forecast4Grid_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for grid data using simplex projection.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbors to use for prediction.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given embedding dimension.
 */
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& source,
                                              const std::vector<std::vector<double>>& target,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<int>& E,
                                              const std::vector<int>& b,
                                              int tau,
                                              int threads);

/*
 * Evaluates prediction performance of different theta parameters for grid data using the S-mapping method.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& source,
                                           const std::vector<std::vector<double>>& target,
                                           const std::vector<int>& lib_indices,
                                           const std::vector<int>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E,
                                           int tau,
                                           int b,
                                           int threads);

/**
 * @brief Evaluate intersection cardinality (IC) for spatial grid embeddings.
 *
 * This function computes the intersection cardinality between the k-nearest neighbors
 * of grid-embedded source and target spatial variables, across a range of embedding dimensions (E)
 * and neighborhood sizes (b). The result is an AUC (Area Under the Curve) score for each (E, b) pair
 * that quantifies the directional similarity or interaction between the spatial fields.
 *
 * The method constructs delay-like embeddings over grid cells using spatial neighborhoods,
 * filters out invalid prediction locations (e.g., with all NaN values), computes nearest neighbors
 * in embedding space, and calculates the cardinality of overlapping neighbors. These overlaps are
 * then evaluated using a CMC-based statistical test (via AUC).
 *
 * Supports both single-threaded and parallel execution using `RcppThread`.
 *
 * @param source 2D spatial variable (grid) used as the source for embedding.
 * @param target 2D spatial variable (grid) used as the target for embedding.
 * @param lib_indices Indices of spatial locations used as the library set (training).
 * @param pred_indices Indices of spatial locations used as the prediction set (evaluation).
 * @param E Vector of spatial embedding dimensions to evaluate (e.g., neighborhood sizes).
 * @param b Vector of neighbor counts (k) used to compute IC.
 * @param tau Spatial embedding spacing (lag). Determines distance between embedding neighbors.
 * @param exclude Number of nearest neighbors to exclude in IC computation.
 * @param threads Maximum number of threads to use.
 * @param parallel_level If > 0, enables parallel evaluation of b for each E.
 *
 * @return A matrix of size (|E| × |b|) × 4 with rows: [E, b, AUC, P-value]
 */
std::vector<std::vector<double>> IC4Grid(const std::vector<std::vector<double>>& source,
                                         const std::vector<std::vector<double>>& target,
                                         const std::vector<size_t>& lib_indices,
                                         const std::vector<size_t>& pred_indices,
                                         const std::vector<int>& E,
                                         const std::vector<int>& b,
                                         int tau,
                                         int exclude,
                                         int threads,
                                         int parallel_level);

#endif // Forecast4Grid_H

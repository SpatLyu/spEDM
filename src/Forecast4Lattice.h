#ifndef Forecast4Lattice_H
#define Forecast4Lattice_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for lattice data using simplex projection.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - nb_vec: A 2D vector of neighbor indices.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbor values to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given combination of E and b.
 */
std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& source,
                                                 const std::vector<double>& target,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<int>& lib_indices,
                                                 const std::vector<int>& pred_indices,
                                                 const std::vector<int>& E,
                                                 const std::vector<int>& b,
                                                 int tau,
                                                 int threads);

/*
 * Evaluates prediction performance of different theta parameters for lattice data using the s-mapping method.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - nb_vec: A 2D vector of neighbor indices.
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
std::vector<std::vector<double>> SMap4Lattice(const std::vector<double>& source,
                                              const std::vector<double>& target,
                                              const std::vector<std::vector<int>>& nb_vec,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<double>& theta,
                                              int E,
                                              int tau,
                                              int b,
                                              int threads);

/**
 * Compute Intersection Cardinality AUC over spatial lattice data.
 *
 * This function computes the causal strength between two lattice-structured time series
 * (`source` and `target`) by evaluating the Intersection Cardinality (IC) curve, and
 * summarizing it using the Area Under the Curve (AUC) metric.
 *
 * For each combination of embedding dimension `E` and neighbor size `b`, the function:
 *  - Generates state-space embeddings based on lattice neighborhood topology.
 *  - Filters out prediction points with missing (NaN) values.
 *  - Computes neighbor structures and evaluates intersection sizes between the mapped
 *    neighbors of `source` and `target`.
 *  - Aggregates the IC curve and estimates the AUC (optionally using significance test).
 *
 * @param source         Time series values of the potential cause variable (flattened lattice vector).
 * @param target         Time series values of the potential effect variable (same shape as `source`).
 * @param nb_vec         Neighborhood topology vector for the lattice structure.
 * @param lib_indices    Indices used for library (training) data.
 * @param pred_indices   Indices used for prediction (testing) data.
 * @param E              Vector of embedding dimensions to try.
 * @param b              Vector of neighbor sizes to try.
 * @param tau            Embedding delay (usually 1 for lattice).
 * @param exclude        Number of nearest neighbors to exclude (e.g., temporal or spatial proximity).
 * @param threads        Number of threads for parallel computation.
 * @param parallel_level Flag indicating whether to use multi-threading (0: serial, 1: parallel).
 *
 * @return A vector of size `E.size() * b.size()`, each element is a vector:
 *         [embedding_dimension, neighbor_size, auc_value, p value].
 *         If inputs are invalid or no prediction point is valid, the AUC value is NaN.
 *
 * @note
 *   - Only AUC and p value are returned in current version. Use other utilities to derive CI.
 *   - Library and prediction indices should be adjusted for 0-based indexing before calling.
 *   - Lattice embedding assumes neighborhood-based spatial structure.
 */
std::vector<std::vector<double>> IC4Lattice(const std::vector<double>& source,
                                            const std::vector<double>& target,
                                            const std::vector<std::vector<int>>& nb_vec,
                                            const std::vector<size_t>& lib_indices,
                                            const std::vector<size_t>& pred_indices,
                                            const std::vector<int>& E,
                                            const std::vector<int>& b,
                                            int tau,
                                            int exclude,
                                            int threads,
                                            int parallel_level);

#endif // Forecast4Lattice_H

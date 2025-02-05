#ifndef SMap_H
#define SMap_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/*
 * Computes the 'S-Maps' prediction.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Spatial cross-section series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
 *   - pred_indices: Vector of T/F values (which states to predict from).
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: A vector<double> containing the predicted target values.
 */
std::vector<double> SMapPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors,
    double theta
);

/*
 * Computes the Rho value using the 'S-Maps' prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Spatial cross-section series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
 *   - pred_indices: Vector of T/F values (which states to predict from).
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: The Pearson correlation coefficient (Rho) between predicted and actual values.
 */
double SMap(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors,
    double theta
);

/*
 * Description: Computes the S-Map prediction and returns a vector containing
 *              Pearson correlation coefficient (PearsonCor), mean absolute error (MAE),
 *              and root mean squared error (RMSE).
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Spatial cross-section series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
 *   - pred_indices: Vector of T/F values (which states to predict from).
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: A vector<double> containing {PearsonCor, MAE, RMSE}.
 */
std::vector<double> SMapBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors,
    double theta
);

#endif // SMap_H

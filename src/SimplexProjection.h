#ifndef SimplexProjection_H
#define SimplexProjection_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

// Description: Computes predictions using the 'simplex projection' method based on state-space reconstruction.
// Parameters:
//   - vectors: Reconstructed state-space (each row represents a separate vector/state).
//   - target: Spatial cross-section series used as the target (should align with vectors).
//   - lib_indices: Vector of T/F values indicating which states to include when searching for neighbors.
//   - pred_indices: Vector of T/F values indicating which states to predict from.
//   - num_neighbors: Number of neighbors to use for simplex projection.
// Returns: A vector<double> containing predicted target values (target_pred).
std::vector<double> SimplexProjectionPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
);

// Description: Computes the Pearson correlation coefficient (rho) using the 'simplex projection' prediction method.
// Parameters:
//   - vectors: Reconstructed state-space (each row represents a separate vector/state).
//   - target: Spatial cross-section series used as the target (should align with vectors).
//   - lib_indices: Vector of T/F values indicating which states to include when searching for neighbors.
//   - pred_indices: Vector of T/F values indicating which states to use for prediction.
//   - num_neighbors: Number of neighbors to use for simplex projection.
// Returns: A double representing the Pearson correlation coefficient (rho) between the predicted and actual target values.
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
);

// Description: Computes the simplex projection and returns a vector containing
//              Pearson correlation coefficient (PearsonCor), mean absolute error (MAE),
//              and root mean squared error (RMSE).
// Parameters:
//   - vectors: Reconstructed state-space (each row is a separate vector/state).
//   - target: Spatial cross-section series to be used as the target (should align with vectors).
//   - lib_indices: Vector of T/F values (which states to include when searching for neighbors).
//   - pred_indices: Vector of T/F values (which states to predict from).
//   - num_neighbors: Number of neighbors to use for simplex projection.
// Returns: A vector<double> containing {PearsonCor, MAE, RMSE}.
std::vector<double> SimplexBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    int num_neighbors
);

#endif // SimplexProjection_H

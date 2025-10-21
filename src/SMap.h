#ifndef SMap_H
#define SMap_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/**
 * @brief Perform S-Mapping prediction using locally weighted linear regression.
 *
 * This function performs prediction based on a reconstructed state-space (spatial-lag embedding).
 * For each prediction index, it:
 *   - Finds the nearest neighbors from the library indices, excluding the current prediction index.
 *   - Computes distance-based weights using the S-mapping weighting parameter (theta).
 *   - Constructs a locally weighted linear regression model using the valid neighbors.
 *   - Predicts the target value using the derived local model.
 *
 * @param vectors        A 2D matrix where each row is a reconstructed state vector (embedding).
 * @param target         A vector of scalar values to predict (e.g., spatial cross sections observations).
 * @param lib_indices    Indices of the vectors used as the library (neighbor candidates).
 * @param pred_indices   Indices of the vectors used for prediction.
 * @param num_neighbors  Number of nearest neighbors to use in local regression. Default is 4.
 * @param theta          Weighting parameter controlling exponential decay of distances. Default is 1.0.
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 * @param dist_average   Whether to average distance by the number of valid vector components. Default is true.
 * @return std::vector<double> Predicted values aligned with the input target vector.
 *         Entries at non-prediction indices or with insufficient valid neighbors are NaN.
 */
std::vector<double> SMapPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/**
 * @brief Perform Composite S-Mapping prediction using multiple reconstructed embeddings.
 *
 * This function extends the standard S-Mapping algorithm to handle multiple
 * reconstructed embeddings (subsets) for each spatial unit. Each subset represents
 * a different reconstructed state-space or embedding dimension group.
 *
 * For each prediction index:
 *   - Computes distances to all library states across each subset independently.
 *   - Averages subset distances (excluding NaN components and subsets with missing data).
 *   - Computes S-map exponential weights based on averaged distances.
 *   - Constructs a locally weighted linear regression model using nearest neighbors.
 *   - Predicts the target value using the fitted local model.
 *
 * Distance calculations exclude NaN components for numerical stability.
 * Supports L1 (Manhattan) or L2 (Euclidean) metrics and optional distance averaging
 * by the number of valid vector components.
 *
 * Supported distance metrics:
 *   dist_metric = 1 → L1 (Manhattan)
 *   dist_metric = 2 → L2 (Euclidean)
 *
 * @param vectors        3D vector of reconstructed embeddings:
 *                       vectors[s][i][j] corresponds to subset s, spatial unit i, and embedding dimension j.
 * @param target         Vector of scalar target values to be predicted.
 * @param lib_indices    Indices of library states (neighbor candidates).
 * @param pred_indices   Indices of prediction targets.
 * @param num_neighbors  Number of nearest neighbors used in local regression (default = 4).
 * @param theta          Weighting parameter controlling exponential decay of distances (default = 1.0).
 * @param dist_metric    Distance metric (1 = L1, 2 = L2). Default = 2.
 * @param dist_average   Whether to average distance by the number of valid vector components (default = true).
 *
 * @return std::vector<double> Predicted values aligned with input target size.
 *         Entries at non-prediction indices or with insufficient neighbors are NaN.
 */
std::vector<double> SMapPrediction(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the Rho value using the 'S-Mapping' prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map. Default is 4.
 *   - theta: Weighting parameter for distances. Default is 1.0.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *
 * Returns: The Pearson correlation coefficient (Rho) between predicted and actual values.
 */
double SMap(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the Rho value using the 'S-Mapping' prediction method (3D version).
 *
 * Each element of vectors is itself a 2D matrix (e.g., multi-component embeddings).
 * The function averages across sub-embeddings before computing distances and predictions.
 */
double SMap(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the S-Mapping prediction and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map. Default is 4.
 *   - theta: Weighting parameter for distances. Default is 1.0.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *
 * Returns: A vector<double> containing {Pearson correlation, MAE, RMSE}.
 */
std::vector<double> SMapBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the S-Mapping prediction and evaluates prediction performance (3D version).
 *
 * Each element of vectors is itself a 2D matrix (e.g., multi-component embeddings).
 * The function averages across sub-embeddings before computing distances and predictions.
 */
std::vector<double> SMapBehavior(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

#endif // SMap_H

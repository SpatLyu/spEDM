#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "NumericUtils.h"
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
) {
  size_t N = target.size();
  std::vector<double> pred(N, std::numeric_limits<double>::quiet_NaN());

  if (num_neighbors <= 0 || lib_indices.empty() || pred_indices.empty()) {
    return pred;
  }

  for (int pred_i : pred_indices) {
    // // Skip if target at prediction index is NaN
    // if (std::isnan(target[pred_i])) {
    //   continue;
    // }

    // Compute distances only for valid vector pairs with valid target values
    std::vector<double> distances;
    std::vector<int> valid_libs;
    distances.reserve(lib_indices.size());
    valid_libs.reserve(lib_indices.size());

    for (int i : lib_indices) {
      // // Only use neighbors with valid target
      // if (std::isnan(target[i])) continue;

      if (i == pred_i) continue; // Skip self-matching

      double sum_sq = 0.0;
      double count = 0.0;
      for (size_t j = 0; j < vectors[pred_i].size(); ++j) {
        double vi = vectors[i][j];
        double vj = vectors[pred_i][j];
        if (!std::isnan(vi) && !std::isnan(vj)) {
          double diff = vi - vj;
          // sum_sq += (dist_metric == 1) ? std::abs(diff) : diff * diff;
          if (dist_metric == 1) {
            sum_sq += std::abs(diff); // L1
          } else {
            sum_sq += diff * diff;    // L2
          }
          count += 1.0;
        }
      }

      if (count > 0) {
        if (dist_metric == 1) {  // L1
          distances.push_back(sum_sq / (dist_average ? count : 1.0));
        } else {                 // L2
          distances.push_back(std::sqrt(sum_sq / (dist_average ? count : 1.0)));
        }
        valid_libs.push_back(i);
      }
    }

    if (distances.empty()) {
      continue; // no usable neighbors
    }

    size_t actual_neighbors = std::min(static_cast<size_t>(num_neighbors), distances.size());

    // Compute mean distance
    double mean_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();

    // Compute weights using exponential kernel
    std::vector<double> weights(distances.size(), 0.0);
    for (size_t i = 0; i < distances.size(); ++i) {
      weights[i] = std::exp(-theta * distances[i] / mean_distance);
    }

    // Select top-k neighbors using partial sort
    std::vector<size_t> neighbor_indices(distances.size());
    std::iota(neighbor_indices.begin(), neighbor_indices.end(), 0);
    std::partial_sort(
      neighbor_indices.begin(),
      neighbor_indices.begin() + actual_neighbors,
      neighbor_indices.end(),
      [&](size_t a, size_t b) {
        return (distances[a] < distances[b]) ||
          (doubleNearlyEqual(distances[a],distances[b]) && a < b);
      }
    );

    // Construct weighted linear system A * coeff = b
    size_t dim = vectors[pred_i].size();
    std::vector<std::vector<double>> A(actual_neighbors, std::vector<double>(dim + 1, 0.0));
    std::vector<double> b(actual_neighbors, 0.0);

    for (size_t i = 0; i < actual_neighbors; ++i) {
      int idx = valid_libs[neighbor_indices[i]];
      double w = weights[neighbor_indices[i]];
      for (size_t j = 0; j < dim; ++j) {
        A[i][j] = vectors[idx][j] * w;
      }
      A[i][dim] = w;  // bias term
      b[i] = target[idx] * w;
    }

    // Solve the system using SVD
    std::vector<std::vector<std::vector<double>>> svd_result = CppSVD(A);
    std::vector<std::vector<double>> U = svd_result[0];
    std::vector<double> S = svd_result[1][0];
    std::vector<std::vector<double>> V = svd_result[2];

    // Invert singular values with tolerance
    double max_s = *std::max_element(S.begin(), S.end());
    std::vector<double> S_inv(S.size(), 0.0);
    for (size_t i = 0; i < S.size(); ++i) {
      if (S[i] >= max_s * 1e-5) {
        S_inv[i] = 1.0 / S[i];
      }
    }

    // Compute regression coefficients: V * S_inv * U^T * b
    std::vector<double> coeff(dim + 1, 0.0);
    for (size_t k = 0; k < V.size(); ++k) {
      double temp = 0.0;
      for (size_t j = 0; j < S_inv.size(); ++j) {
        for (size_t i = 0; i < U.size(); ++i) {
          temp += V[k][j] * S_inv[j] * U[i][j] * b[i];
        }
      }
      coeff[k] = temp;
    }

    // Compute prediction: dot(coeff, input) + bias
    double prediction = 0.0;
    for (size_t i = 0; i < dim; ++i) {
      prediction += coeff[i] * vectors[pred_i][i];
    }
    prediction += coeff[dim];

    pred[pred_i] = prediction;
  }

  return pred;
}

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
) {
  size_t N = target.size();
  std::vector<double> pred(N, std::numeric_limits<double>::quiet_NaN());

  if (num_neighbors <= 0 || lib_indices.empty() || pred_indices.empty() || vectors.empty()) {
    return pred;
  }

  size_t num_subsets = vectors.size();

  for (int pred_i : pred_indices) {
    if (pred_i < 0 || static_cast<size_t>(pred_i) >= N) continue;

    // Compute averaged distances across all subsets
    std::vector<double> distances;
    std::vector<int> valid_libs;
    distances.reserve(lib_indices.size());
    valid_libs.reserve(lib_indices.size());

    for (int i : lib_indices) {
      if (i == pred_i || i < 0 || static_cast<size_t>(i) >= N) continue;

      std::vector<double> subset_distances;
      subset_distances.reserve(num_subsets);

      // Compute distance within each subset separately
      for (size_t s = 0; s < num_subsets; ++s) {
        if (i >= static_cast<int>(vectors[s].size()) || pred_i >= static_cast<int>(vectors[s].size()))
          continue;

        const auto& vec_i = vectors[s][i];
        const auto& vec_p = vectors[s][pred_i];

        double sum_sq = 0.0;
        double count = 0.0;

        for (size_t j = 0; j < vec_p.size(); ++j) {
          if (j < vec_i.size() && !std::isnan(vec_i[j]) && !std::isnan(vec_p[j])) {
            double diff = vec_i[j] - vec_p[j];
            if (dist_metric == 1) {
              sum_sq += std::abs(diff); // L1
            } else {
              sum_sq += diff * diff;    // L2
            }
            count += 1.0;
          }
        }

        if (count > 0) {
          if (dist_metric == 1) {
            subset_distances.push_back(sum_sq / (dist_average ? count : 1.0));
          } else {
            subset_distances.push_back(std::sqrt(sum_sq / (dist_average ? count : 1.0)));
          }
        }
      }

      // Average distances across subsets (exclude NaNs)
      if (!subset_distances.empty()) {
        double avg_dist = std::accumulate(subset_distances.begin(), subset_distances.end(), 0.0) /
          subset_distances.size();
        distances.push_back(avg_dist);
        valid_libs.push_back(i);
      }
    }

    if (distances.empty()) {
      continue; // No valid distances found
    }

    size_t actual_neighbors = std::min(static_cast<size_t>(num_neighbors), distances.size());

    // Compute mean distance for weighting
    double mean_distance = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();

    // Compute S-map exponential weights
    std::vector<double> weights(distances.size(), 0.0);
    for (size_t i = 0; i < distances.size(); ++i) {
      weights[i] = std::exp(-theta * distances[i] / mean_distance);
    }

    // Select top-k nearest neighbors
    std::vector<size_t> neighbor_indices(distances.size());
    std::iota(neighbor_indices.begin(), neighbor_indices.end(), 0);
    std::partial_sort(
      neighbor_indices.begin(),
      neighbor_indices.begin() + actual_neighbors,
      neighbor_indices.end(),
      [&](size_t a, size_t b) {
        return (distances[a] < distances[b]) ||
          (doubleNearlyEqual(distances[a],distances[b]) && a < b);
      }
    );

    // Build weighted linear regression system
    // Concatenate embeddings from all subsets for each vector
    size_t dim = 0;
    for (const auto& subset : vectors) {
      if (!subset.empty()) dim += subset[0].size();
    }

    std::vector<std::vector<double>> A(actual_neighbors, std::vector<double>(dim + 1, 0.0));
    std::vector<double> b(actual_neighbors, 0.0);

    for (size_t ni = 0; ni < actual_neighbors; ++ni) {
      int idx = valid_libs[neighbor_indices[ni]];
      double w = weights[neighbor_indices[ni]];

      // Flatten all subsets into a single combined vector for regression
      size_t pos = 0;
      for (size_t s = 0; s < num_subsets; ++s) {
        if (idx >= static_cast<int>(vectors[s].size())) continue;
        const auto& vec = vectors[s][idx];
        for (double val : vec) {
          A[ni][pos++] = val * w;
        }
      }

      A[ni][dim] = w;           // bias term
      b[ni] = target[idx] * w;  // weighted target
    }

    // Solve via SVD: V * S⁻¹ * Uᵀ * b
    std::vector<std::vector<std::vector<double>>> svd_result = CppSVD(A);
    std::vector<std::vector<double>> U = svd_result[0];
    std::vector<double> S = svd_result[1][0];
    std::vector<std::vector<double>> V = svd_result[2];

    double max_s = *std::max_element(S.begin(), S.end());
    std::vector<double> S_inv(S.size(), 0.0);
    for (size_t i = 0; i < S.size(); ++i) {
      if (S[i] >= max_s * 1e-5) {
        S_inv[i] = 1.0 / S[i];
      }
    }

    std::vector<double> coeff(dim + 1, 0.0);
    for (size_t k = 0; k < V.size(); ++k) {
      double temp = 0.0;
      for (size_t j = 0; j < S_inv.size(); ++j) {
        for (size_t i = 0; i < U.size(); ++i) {
          temp += V[k][j] * S_inv[j] * U[i][j] * b[i];
        }
      }
      coeff[k] = temp;
    }

    // Compute prediction using flattened input vector
    double prediction = 0.0;
    size_t pos = 0;
    for (size_t s = 0; s < num_subsets; ++s) {
      if (pred_i >= static_cast<int>(vectors[s].size())) continue;
      for (double val : vectors[s][pred_i]) {
        prediction += coeff[pos++] * val;
      }
    }
    prediction += coeff[dim]; // bias

    pred[pred_i] = prediction;
  }

  return pred;
}

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
) {
  double rho = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta, dist_metric, dist_average);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    rho = PearsonCor(target_pred, target, true);
  }
  return rho;
}

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
) {
  double rho = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta, dist_metric, dist_average);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    rho = PearsonCor(target_pred, target, true);
  }
  return rho;
}

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
) {
  // Initialize PearsonCor, MAE, and RMSE
  double pearson = std::numeric_limits<double>::quiet_NaN();
  double mae = std::numeric_limits<double>::quiet_NaN();
  double rmse = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta, dist_metric, dist_average);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    // Compute PearsonCor, MAE, and RMSE
    pearson = PearsonCor(target_pred, target, true);
    mae = CppMAE(target_pred, target, true);
    rmse = CppRMSE(target_pred, target, true);
  }

  // Return the three metrics as a vector
  return {pearson, mae, rmse};
}

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
) {
  // Initialize PearsonCor, MAE, and RMSE
  double pearson = std::numeric_limits<double>::quiet_NaN();
  double mae = std::numeric_limits<double>::quiet_NaN();
  double rmse = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta, dist_metric, dist_average);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    // Compute PearsonCor, MAE, and RMSE
    pearson = PearsonCor(target_pred, target, true);
    mae = CppMAE(target_pred, target, true);
    rmse = CppRMSE(target_pred, target, true);
  }

  // Return the three metrics as a vector
  return {pearson, mae, rmse};
}

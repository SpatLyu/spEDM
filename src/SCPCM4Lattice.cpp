#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "spEDMDataStruct.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the partial correlation between the target variable and its simplex projection,
 *        incorporating control variables using a lattice-based embedding approach.
 *
 * @param vectors: Reconstructed state-space, where each row represents a separate state vector.
 * @param target: Spatial cross-section series to be used as the target, aligned with 'vectors'.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param nb_vec: Neighbor indices for each spatial unit.
 * @param lib_indices: Boolean vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Boolean vector indicating which states to use for predictions.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for simplex projection.
 * @param cumulate: Flag indicating whether to cumulatively incorporate control variables.
 *
 * @return A std::vector<double> containing:
 *         - rho[0]: Pearson correlation between the target and its simplex projection.
 *         - rho[1]: Partial correlation controlling for the influence of the control variables.
 */
std::vector<double> PartialSimplex4Lattice(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<std::vector<int>>& nb_vec,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate) {
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      } else {
        temp_pred = SimplexProjectionPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i]);
      }
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls]);
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],taus[i]);
      temp_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1]);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCor(target,target_pred,con_pred,true,false);
  }
  return rho;
}

/**
 * @brief Computes the partial correlation between a spatial cross-section series and its prediction
 *        using the S-Map method, incorporating control variables.
 *
 * This function performs state-space reconstruction and S-Map prediction while accounting for
 * control variables in a lattice-based spatial setting. The process can be either cumulative or
 * independent in terms of incorporating control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Spatial cross-section series used as the prediction target.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param nb_vec: Neighbor indices vector specifying spatial unit neighbors.
 * @param lib_indices: Boolean vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Boolean vector indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for S-Map prediction.
 * @param theta: Weighting parameter for distances in S-Map.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSMap4Lattice(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<std::vector<int>>& nb_vec,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    double theta,
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate){
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      } else {
        temp_pred = SMapPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i], theta);
      }
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls], theta);
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],taus[i]);
      temp_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1], theta);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCor(target,target_pred,con_pred,true,false);
  }

  return rho;
}

/*
 * Perform SCPCM on a single library and prediction set for lattice data.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: Spatial cross-section series used as the target (should align with x_vectors).
 *   - controls: Cross-sectional data of control variables (stored by row).
 *   - nb_vec: Neighbor indices vector of the spatial units.
 *   - lib_indices: A boolean vector indicating which states to include when searching for neighbors.
 *   - lib_size: Size of the library used for cross mapping.
 *   - max_lib_size: Maximum size of the library.
 *   - possible_lib_indices: Indices of possible library states.
 *   - pred_indices: A boolean vector indicating which states to use for prediction.
 *   - conEs: Number of dimensions for attractor reconstruction with control variables.
 *   - taus: Spatial lag step for constructing lagged state-space vectors with control variables.
 *   - b: A vector specifying the numbers of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *   - cumulate: Whether to accumulate partial correlations.
 *
 * Returns:
 *   A vector of PartialCorRes objects, where each contains:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho).
 *   - A double representing the Partial correlation coefficient (pratial rho).
 */
std::vector<PartialCorRes> SCPCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                       // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,        // Neighbor indices vector of the spatial units
    const std::vector<bool>& lib_indices,               // Vector of T/F values (which states to include when searching for neighbors)
    int lib_size,                                       // Size of the library
    int max_lib_size,                                   // Maximum size of the library
    const std::vector<int>& possible_lib_indices,       // Indices of possible library states
    const std::vector<bool>& pred_indices,              // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                      // Number of dimensions for the attractor reconstruction with control variables
    const std::vector<int>& taus,                       // Spatial lag step for constructing lagged state-space vectors with control variables
    const std::vector<int>& b,                          // Numbers of neighbors to use for simplex projection
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    bool cumulate                                       // Whether to cumulate the partial correlations
) {
  int n = x_vectors.size();
  std::vector<PartialCorRes> x_xmap_y;
  std::vector<double> rho;

  if (lib_size == max_lib_size) { // No possible library variation if using all vectors
    std::vector<bool> lib_indices(n, false);
    for (int idx : possible_lib_indices) {
      lib_indices[idx] = true;
    }

    // Run partial cross map and store results
    if (simplex) {
      rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, taus, b, cumulate);
    } else {
      rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, taus, b, theta, cumulate);
    }
    x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
  } else {
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<bool> lib_indices(n, false);
      // Setup changing library
      if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
        for (int i = start_lib; i < max_lib_size; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          lib_indices[possible_lib_indices[i]] = true;
        }
      }

      // Run partial cross map and store results
      if (simplex) {
        rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, taus, b, cumulate);
      } else {
        rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, taus, b, theta, cumulate);
      }
      x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    }
  }

  return x_xmap_y;
}

/**
 * Performs SCPCM on a spatial lattice dataset.
 *
 * Parameters:
 * - x: Spatial cross-section series used as the predictor variable (**cross mapping from**).
 * - y: Spatial cross-section series used as the target variable (**cross mapping to**).
 * - controls: Cross-sectional data of control variables (**stored by row**).
 * - nb_vec: A nested vector containing neighborhood information for lattice data.
 * - lib_sizes: A vector specifying different library sizes for SCPCM analysis.
 * - lib: A vector specifying the library indices (1-based in R, converted to 0-based in C++).
 * - pred: A vector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 * - Es: A vector specifying the embedding dimensions for attractor reconstruction using x and control variables.
 * - taus: A vector specifying the spatial lag steps for constructing lagged state-space vectors using x and control variables.
 * - b: A vector specifying the numbers of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - cumulate: Boolean flag indicating whether to cumulate partial correlations.
 * - progressbar: Boolean flag indicating whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean pearson cross-mapping correlation.
 *      - The statistical significance of the pearson correlation.
 *      - The lower bound of the pearson correlation confidence interval.
 *      - The upper bound of the pearson correlation confidence interval.
 *      - The mean partial cross-mapping partial correlation.
 *      - The statistical significance of the partial correlation.
 *      - The lower bound of the partial correlation confidence interval.
 *      - The upper bound of the partial correlation confidence interval.
 */
std::vector<std::vector<double>> SCPCM4Lattice(
    const std::vector<double>& x,                       // Spatial cross-section series to cross map from
    const std::vector<double>& y,                       // Spatial cross-section series to cross map to
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,        // Neighbor indices vector of the spatial units
    const std::vector<int>& lib_sizes,                  // Vector of library sizes to use
    const std::vector<int>& lib,                        // Vector specifying the library indices
    const std::vector<int>& pred,                       // Vector specifying the prediction indices
    const std::vector<int>& Es,                         // Number of dimensions for the attractor reconstruction with the x and control variables
    const std::vector<int>& taus,                       // Spatial lag step for constructing lagged state-space vectors with the x and control variables
    const std::vector<int>& b,                          // Numbers of nearest neighbors to use for prediction
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                        // Number of threads used from the global pool
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool progressbar                                    // Whether to print the progress bar
) {
  // If b is not provided correctly, default it to E + 2
  std::vector<int> bs = b;
  for (size_t i = 0; i < bs.size(); ++i){
    if (bs[i] <= 0) {
      bs[i] = Es[i] + 2;
    }
  }

  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  int taux = taus[0];
  std::vector<int> contaus = taus;
  contaus.erase(contaus.begin());

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(threads);
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  std::vector<std::vector<double>> x_vectors = GenLatticeEmbeddings(x,nb_vec,Ex,taux);
  int n = x_vectors.size();

  int n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = controls.size();
  }

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> lib_indices(n, false);
  std::vector<bool> pred_indices(n, false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int libsize_int = lib.size();
  for (int i = 0; i < libsize_int; ++i) {
    // // Do not strictly exclude spatial units with embedded state-space vectors containing NaN values from participating in cross mapping.
    // if (!checkOneDimVectorHasNaN(x_vectors[lib[i] - 1])){
    //   lib_indices[lib[i] - 1] = true;
    // }
    lib_indices[lib[i] - 1] = true; // Convert to 0-based index
  }
  int predsize_int = pred.size();
  for (int i = 0; i < predsize_int; ++i) {
    // // Do not strictly exclude spatial units with embedded state-space vectors containing NaN values from participating in cross mapping.
    // if (!checkOneDimVectorHasNaN(x_vectors[pred[i] - 1])){
    //   pred_indices[pred[i] - 1] = true;
    // }
    pred_indices[pred[i] - 1] = true; // Convert to 0-based index
  }

  // /* Aligned with the previous implementation,
  //  * now deprecated in the source C++ code.
  //  *  ----- Wenbo Lv, written on 2025.02.10
  //  */
  // for (int i = 0, i < (Ex - 1) * tau, ++i){
  //   lib_indices[lib[i] - 1] = false;
  //   pred_indices[pred[i] - 1] = false;
  // }

  // /* Do not uncomment those codes;
  //  * it's the previous implementation using `std::vector<std::pair<int, int>>`  input for lib and pred,
  //  * kept for reference. ----- Wenbo Lv, written on 2025.02.09
  //  */
  // // Setup pred_indices
  // std::vector<bool> pred_indices(n, false);
  // for (const auto& p : pred) {
  //   int row_start = p.first + (Ex - 1) * tau;
  //   int row_end = p.second;
  //   if (row_end > row_start && row_start >= 0 && row_end < n) {
  //     std::fill(pred_indices.begin() + row_start, pred_indices.begin() + row_end + 1, true);
  //   }
  // }
  //
  // // Setup lib_indices
  // std::vector<bool> lib_indices(n, false);
  // for (const auto& l : lib) {
  //   int row_start = l.first + (Ex - 1) * tau;
  //   int row_end = l.second;
  //   if (row_end > row_start && row_start >= 0 && row_end < n) {
  //     std::fill(lib_indices.begin() + row_start, lib_indices.begin() + row_end + 1, true);
  //   }
  // }

  int max_lib_size = std::accumulate(lib_indices.begin(), lib_indices.end(), 0); // Maximum lib size
  std::vector<int> possible_lib_indices;
  for (int i = 0; i < n; ++i) {
    if (lib_indices[i]) {
      possible_lib_indices.push_back(i);
    }
  }

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is Ex + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, Ex + 2); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Initialize the result container
  std::vector<PartialCorRes> x_xmap_y;

  // Sequential version of the for loop
  // for (int lib_size : unique_lib_sizes) {
  //   RcppThread::Rcout << "lib_size: " << lib_size << "\n";
  //   auto results = SCPCMSingle4Lattice(
  //     x_vectors,
  //     y,
  //     controls,
  //     nb_vec,
  //     lib_indices,
  //     lib_size,
  //     max_lib_size,
  //     possible_lib_indices,
  //     pred_indices,
  //     conEs,
  //     contaus,
  //     bs,
  //     simplex,
  //     theta,
  //     cumulate
  //   );
  //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
  // }

  // Perform the operations using RcppThread
  if (progressbar) {
    RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = SCPCMSingle4Lattice(
        x_vectors,
        y,
        controls,
        nb_vec,
        lib_indices,
        lib_size,
        max_lib_size,
        possible_lib_indices,
        pred_indices,
        conEs,
        contaus,
        bs,
        simplex,
        theta,
        cumulate
      );
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = SCPCMSingle4Lattice(
        x_vectors,
        y,
        controls,
        nb_vec,
        lib_indices,
        lib_size,
        max_lib_size,
        possible_lib_indices,
        pred_indices,
        conEs,
        contaus,
        bs,
        simplex,
        theta,
        cumulate
      );
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    }, threads_sizet);
  }

  // Group by the first int and store second and third values as pairs
  std::map<int, std::vector<std::pair<double, double>>> grouped_results;

  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].emplace_back(result.second, result.third);
  }

  std::vector<std::vector<double>> final_results;

  // Compute the mean of second and third values for each group
  for (const auto& group : grouped_results) {
    std::vector<double> second_values, third_values;

    for (const auto& val : group.second) {
      second_values.push_back(val.first);
      third_values.push_back(val.second);
    }

    double mean_second = CppMean(second_values, true);
    double mean_third = CppMean(third_values, true);

    final_results.push_back({static_cast<double>(group.first), mean_second, mean_third});
  }

  // Compute significance and confidence intervals for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho_second = final_results[i][1];
    double rho_third = final_results[i][2];

    // Compute significance and confidence interval for second value
    double significance_second = CppCorSignificance(rho_second, n);
    std::vector<double> confidence_interval_second = CppCorConfidence(rho_second, n);

    // Compute significance and confidence interval for third value
    double significance_third = CppCorSignificance(rho_third, n, n_confounds);
    std::vector<double> confidence_interval_third = CppCorConfidence(rho_third, n, n_confounds);

    // Append computed statistical values to the result
    final_results[i].push_back(significance_second);
    final_results[i].push_back(confidence_interval_second[0]);
    final_results[i].push_back(confidence_interval_second[1]);

    final_results[i].push_back(significance_third);
    final_results[i].push_back(confidence_interval_third[0]);
    final_results[i].push_back(confidence_interval_third[1]);
  }

  return final_results;
}

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
 * @param num_neighbors: Number of neighbors to use for simplex projection.
 * @param cumulate: Flag indicating whether to cumulatively incorporate control variables.
 * @param includeself: Flag indicating whether to include the current state in the embedding vector.
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
    int num_neighbors,
    bool cumulate,
    bool includeself
){
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate) {
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors);
      } else {
        temp_pred = SimplexProjectionPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors);
      }
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],includeself);
    }

    std::vector<double> con_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors);
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors);
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],includeself);
      temp_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

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
 * @param num_neighbors: Number of neighbors to use for S-Map prediction.
 * @param theta: Weighting parameter for distances in S-Map.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @param includeself: Boolean flag to indicate whether to include the current state in the embedding.
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
    int num_neighbors,
    double theta,
    bool cumulate,
    bool includeself
){
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate){
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      } else {
        temp_pred = SMapPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      }
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],includeself);
    }

    std::vector<double> con_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors, theta);
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      temp_embedding = GenLatticeEmbeddings(temp_pred,nb_vec,conEs[i],includeself);
      temp_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors, theta);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCor(target,target_pred,con_pred,true,false);
  }

  return rho;
}

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
    int E,                                              // Number of dimensions for the attractor reconstruction
    int b,                                              // Number of neighbors to use for simplex projection
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool includeself                                    // Whether to include the current state when constructing the embedding vector
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
      rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, cumulate, includeself);
    } else {
      rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, theta, cumulate, includeself);
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
        rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, cumulate, includeself);
      } else {
        rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, theta, cumulate, includeself);
      }
      x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    }
  }

  return x_xmap_y;
}

std::vector<std::vector<double>> SCPCM4Lattice(
    const std::vector<double>& x,                       // Spatial cross-section series to cross map from
    const std::vector<double>& y,                       // Spatial cross-section series to cross map to
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,        // Neighbor indices vector of the spatial units
    const std::vector<int>& lib_sizes,                  // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& lib,        // Matrix (n x 2) using n sequences of data to construct libraries
    const std::vector<std::pair<int, int>>& pred,       // Matrix (n x 2) using n sequences of data to predict from
    const std::vector<int>& Es,                         // Number of dimensions for the attractor reconstruction with the x and control variables
    int tau,                                            // Spatial lag for the lagged-vector construction
    int b,                                              // Number of nearest neighbors to use for prediction
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                        // Number of threads used from the global pool
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool includeself,                                   // Whether to include the current state when constructing the embedding vector
    bool progressbar = true                             // Whether to print the progress bar
) {
  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  // If b is not provided correctly, default it to Ex + 2
  if (b <= 0) {
    b = Ex + 2;
  }

  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  std::vector<std::vector<double>> x_vectors = GenLatticeEmbeddings(x,nb_vec,Ex,includeself);
  int n = x_vectors.size();

  int n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = controls.size();
  }

  // Setup pred_indices
  std::vector<bool> pred_indices(n, false);
  for (const auto& p : pred) {
    int row_start = p.first + (Ex - 1) * tau;
    int row_end = p.second;
    if (row_end > row_start && row_start >= 0 && row_end < n) {
      std::fill(pred_indices.begin() + row_start, pred_indices.begin() + row_end + 1, true);
    }
  }

  // Setup lib_indices
  std::vector<bool> lib_indices(n, false);
  for (const auto& l : lib) {
    int row_start = l.first + (Ex - 1) * tau;
    int row_end = l.second;
    if (row_end > row_start && row_start >= 0 && row_end < n) {
      std::fill(lib_indices.begin() + row_start, lib_indices.begin() + row_end + 1, true);
    }
  }

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
  //     Ex,
  //     b,
  //     simplex,
  //     theta,
  //     cumulate,
  //     includeself
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
        Ex,
        b,
        simplex,
        theta,
        cumulate,
        includeself
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
        Ex,
        b,
        simplex,
        theta,
        cumulate,
        includeself
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

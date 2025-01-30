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
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

double PartialSimplex4Lattice(
    const std::vector<std::vector<double>>& vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls, // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,      // Neighbor indices vector of the spatial units
    const std::vector<bool>& lib_indices,             // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,            // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                    // Number of dimensions for the attractor reconstruction with control variables
    int num_neighbors,                                // Number of neighbors to use for simplex projection
    bool includeself                                  // Whether to include the current state when constructing the embedding vector
){
  int n_controls = controls.size();
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

  double rho = PartialCorTrivar(target,target_pred,con_pred,true,false);
  return rho;
}

double PartialSMap4Lattice(
    const std::vector<std::vector<double>>& vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls, // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,      // Neighbor indices vector of the spatial units
    const std::vector<bool>& lib_indices,             // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,            // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                    // Number of dimensions for the attractor reconstruction with control variables
    int num_neighbors,                                // Number of neighbors to use for simplex projection
    double theta,                                     // Weighting parameter for distances
    bool includeself                                  // Whether to include the current state when constructing the embedding vector
){
  int n_controls = controls.size();
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

  double rho = PartialCorTrivar(target,target_pred,con_pred,true,false);
  return rho;
}

std::vector<std::pair<int, double>> SCPCMSingle4Lattice(
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
    bool includeself                                    // Whether to include the current state when constructing the embedding vector
) {
  int n = x_vectors.size();
  std::vector<std::pair<int, double>> x_xmap_y;
  double rho;

  if (lib_size == max_lib_size) { // No possible library variation if using all vectors
    std::vector<bool> lib_indices(n, false);
    for (int idx : possible_lib_indices) {
      lib_indices[idx] = true;
    }

    // Run partial cross map and store results
    if (simplex) {
      rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, includeself);
    } else {
      rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, theta, includeself);
    }
    x_xmap_y.emplace_back(lib_size, rho);
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
        rho = PartialSimplex4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, includeself);
      } else {
        rho = PartialSMap4Lattice(x_vectors, y, controls, nb_vec, lib_indices, pred_indices, conEs, b, theta, includeself);
      }
      x_xmap_y.emplace_back(lib_size, rho);
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
  std::vector<std::pair<int, double>> x_xmap_y;

  // // Sequential version of the for loop
  // for (int lib_size : unique_lib_sizes) {
  //   auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b);
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
        includeself
      );
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    }, threads_sizet);
  }

  // Group by the first int and compute the mean
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  std::vector<std::vector<double>> final_results;
  for (const auto& group : grouped_results) {
    double mean_value = CppMean(group.second, true);
    final_results.push_back({static_cast<double>(group.first), mean_value});
  }

  // Calculate significance and confidence interval for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho = final_results[i][1];
    double significance = CppSignificance(rho, n);
    std::vector<double> confidence_interval = CppConfidence(rho, n);

    final_results[i].push_back(significance);
    final_results[i].push_back(confidence_interval[0]);
    final_results[i].push_back(confidence_interval[1]);
  }

  return final_results;
}

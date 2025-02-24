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

/*
 * Perform GCCM on a single lib and pred for lattice data.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: Spatial cross-section series used as the target (should align with x_vectors).
 *   - lib_indices: A boolean vector indicating which states to include when searching for neighbors.
 *   - lib_size: Size of the library used for cross mapping.
 *   - max_lib_size: Maximum size of the library.
 *   - possible_lib_indices: Indices of possible library states.
 *   - pred_indices: A boolean vector indicating which states to use for prediction.
 *   - b: Number of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *
 * Returns:
 *   A vector of pairs, where each pair consists of:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho) between predicted and actual values.
 */
std::vector<std::pair<int, double>> GCCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,
    const std::vector<double>& y,
    const std::vector<bool>& lib_indices,
    int lib_size,
    int max_lib_size,
    const std::vector<int>& possible_lib_indices,
    const std::vector<bool>& pred_indices,
    int b,
    bool simplex,
    double theta
) {
  int n = x_vectors.size();
  std::vector<std::pair<int, double>> x_xmap_y;
  double rho;

  if (lib_size == max_lib_size) { // No possible library variation if using all vectors
    std::vector<bool> lib_indices(n, false);
    for (int idx : possible_lib_indices) {
      lib_indices[idx] = true;
    }

    // Run cross map and store results
    if (simplex) {
      rho = SimplexProjection(x_vectors, y, lib_indices, pred_indices, b);
    } else {
      rho = SMap(x_vectors, y, lib_indices, pred_indices, b, theta);
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

      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(x_vectors, y, lib_indices, pred_indices, b);
      } else {
        rho = SMap(x_vectors, y, lib_indices, pred_indices, b, theta);
      }
      x_xmap_y.emplace_back(lib_size, rho);
    }
  }

  return x_xmap_y;
}

/**
 * Performs GCCM on a spatial lattice data.
 *
 * Parameters:
 * - x: Spatial cross-section series used as the predict variable (**cross mapping from**).
 * - y: Spatial cross-section series used as the target variable (**cross mapping to**).
 * - nb_vec: A nested vector containing neighborhood information for lattice data.
 * - lib_sizes: A vector specifying different library sizes for GCCM analysis.
 * - lib: A vector specifying the library indices (1-based in R, converted to 0-based in C++).
 * - pred: A vector specifying the prediction indices (1-based in R, converted to 0-based in C++).
 * - E: Embedding dimension for attractor reconstruction.
 * - tau: the step of spatial lags for prediction.
 * - b: Number of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - progressbar: Boolean flag to indicate whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean cross-mapping correlation.
 *      - The statistical significance of the correlation.
 *      - The lower bound of the confidence interval.
 *      - The upper bound of the confidence interval.
 */
std::vector<std::vector<double>> GCCM4Lattice(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::vector<int>>& nb_vec,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    bool progressbar
) {
  // If b is not provided correctly, default it to E + 2
  if (b <= 0) {
    b = E + 2;
  }

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(threads);
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  std::vector<std::vector<double>> x_vectors = GenLatticeEmbeddings(x, nb_vec, E, tau);

  int n = x_vectors.size();

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
  // for (int i = 0, i < (E - 1) * tau, ++i){
  //   lib_indices[lib[i] - 1] = false;
  //   pred_indices[pred[i] - 1] = false;
  // }

  // /* Do not uncomment those codes;
  //  * it's the previous implementation using `std::vector<std::pair<int, int>>`  input for lib and pred,
  //  * kept for reference. ----- Wenbo Lv, written on 2025.02.09
  //  */
  // // Setup pred_indices
  // std::vector<bool> pred_indices(n, false);
  //
  // for (const auto& p : pred) {
  //   int row_start = p.first + (E - 1) * tau;
  //   int row_end = p.second;
  //   if (row_end > row_start && row_start >= 0 && row_end < n) {
  //     std::fill(pred_indices.begin() + row_start, pred_indices.begin() + row_end + 1, true);
  //   }
  // }
  //
  // // Setup lib_indices
  // std::vector<bool> lib_indices(n, false);
  // for (const auto& l : lib) {
  //   int row_start = l.first + (E - 1) * tau;
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

  // Ensure the minimum value in unique_lib_sizes is E + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, E + 2); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // // Sequential version of the for loop
  // for (int lib_size : unique_lib_sizes) {
  //   auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b, simplex, theta);
  //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
  // }

  // Perform the operations using RcppThread
  if (progressbar) {
    RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b, simplex, theta);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = GCCMSingle4Lattice(x_vectors, y, lib_indices, lib_size, max_lib_size, possible_lib_indices, pred_indices, b, simplex, theta);
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
    double significance = CppCorSignificance(rho, n);
    std::vector<double> confidence_interval = CppCorConfidence(rho, n);

    final_results[i].push_back(significance);
    final_results[i].push_back(confidence_interval[0]);
    final_results[i].push_back(confidence_interval[1]);
  }

  return final_results;
}

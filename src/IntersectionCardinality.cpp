#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include "CppStats.h"
#include "spEDMDataStruct.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * Computes intersection-based mapping ratio sequences between two neighbor graphs
 * for use in Cross Mapping Cardinality (CMC) or similar causal inference frameworks.
 *
 * Parameters:
 *   neighborsX     - Precomputed sorted neighbor indices for embedding X
 *   neighborsY     - Precomputed sorted neighbor indices for embedding Y
 *   lib_size       - Size of the moving library used in mapping
 *   lib_indices    - Global indices from which to draw the sliding libraries
 *   pred_indices   - Indices at which to perform prediction (evaluation points)
 *   num_neighbors  - Number of neighbors used for mapping (after exclusion)
 *   n_excluded     - Number of nearest neighbors to exclude from the front
 *   threads        - Number of parallel threads for computation
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of IntersectionRes structures, each containing the average intersection
 *   ratio sequence (IC curve) for a different starting point of the moving library.
 *   If lib_size == lib_indices.size(), returns a single result using full library.
 *
 * Notes:
 *   - Neighbor lists must use std::numeric_limits<size_t>::max() to indicate invalid entries.
 *   - This function assumes that the neighbor vectors are sorted by ascending distance.
 *   - Use in combination with AUC computation to assess causal strength.
 */
std::vector<IntersectionRes> IntersectionCardinalitySingle(
    const std::vector<std::vector<size_t>>& neighborsX,
    const std::vector<std::vector<size_t>>& neighborsY,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    size_t num_neighbors,
    size_t n_excluded,
    size_t threads,
    int parallel_level
){
  int max_lib_size = lib_indices.size();
  const size_t max_r = num_neighbors + n_excluded; // Total number of neighbors = actual used + excluded ones

  auto ICSingle = [&](const std::vector<int>& lib) {
    // Store mapping ratio curves for each prediction point
    std::vector<std::vector<double>> ratio_curves(pred_indices.size(), std::vector<double>(num_neighbors, std::numeric_limits<double>::quiet_NaN()));

    // Precompute library set for fast lookup
    std::unordered_set<size_t> lib_set(lib.begin(), lib.end());

    if (parallel_level == 0){
      // Perform the operations using RcppThread
      RcppThread::parallelFor(0, pred_indices.size(), [&](size_t i) {
        const int idx = pred_indices[i];

        if ((neighborsX[idx][0] != std::numeric_limits<size_t>::max()) &&
            (neighborsY[idx][0] != std::numeric_limits<size_t>::max())){
          if ((neighborsX[idx].size() >= max_r) && (neighborsY[idx].size() >= max_r)){
            // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_x;
            for (size_t iidx = 0; iidx < neighborsX[idx].size(); ++iidx){
              if(lib_set.count(neighborsX[idx][iidx]) > 0){
                neighbors_x.push_back(neighborsX[idx][iidx]);
                if (neighbors_x.size() > max_r) break;
              }
            }
            if (neighbors_x.size() > n_excluded) {
              neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
            }

            // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_y;
            for (size_t iidx = 0; iidx < neighborsY[idx].size(); ++iidx){
              if(lib_set.count(neighborsY[idx][iidx])){
                neighbors_y.push_back(neighborsY[idx][iidx]);
                if (neighbors_y.size() > max_r) break;
              }
            }
            if (neighbors_y.size() > n_excluded) {
              neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
            }

            // Precompute y-neighbors set for fast lookup
            std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

            // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
            std::unordered_map<size_t, std::vector<size_t>> mapped_neighbors;

            for (size_t nx : neighbors_x) {
              if (neighborsY[nx][0] != std::numeric_limits<size_t>::max()) {
                for (size_t iidx = 0; iidx < neighborsY[nx].size(); ++iidx) {
                  if (lib_set.count(neighborsY[nx][iidx])) {
                    mapped_neighbors[nx].push_back(neighborsY[nx][iidx]);
                    if (mapped_neighbors[nx].size() >= num_neighbors) break;
                  }
                }
              }
            }

            // Compute intersection ratio between mapped x-neighbors and original y-neighbors
            for (size_t ki = 0; ki < num_neighbors; ++ki) {
              size_t count = 0;

              for (size_t nx : neighbors_x) {
                auto it = mapped_neighbors.find(nx);
                if (it != mapped_neighbors.end() && ki < it->second.size()) {
                  const auto& yn = it->second;

                  // Check if any of the first (ki+1) mapped neighbors exist in y's original neighbors
                  for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
                    if (y_neighbors_set.count(yn[pos])) {
                      ++count;
                      break;  // Count each x-neighbor only once if it intersects
                    }
                  }
                }
              }

              if (neighbors_x.size() > 0) {
                ratio_curves[i][ki] = static_cast<double>(count) / static_cast<double>(neighbors_x.size());
              }
            }
          }
        }
      }, threads);
    } else {
      // Perform the operations one by one
      for (size_t i = 0; i < pred_indices.size(); ++i){
        const int idx = pred_indices[i];

        if ((neighborsX[idx][0] != std::numeric_limits<size_t>::max()) &&
            (neighborsY[idx][0] != std::numeric_limits<size_t>::max())){
          if ((neighborsX[idx].size() >= max_r) && (neighborsY[idx].size() >= max_r)){
            // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_x;
            for (size_t iidx = 0; iidx < neighborsX[idx].size(); ++iidx){
              if(lib_set.count(neighborsX[idx][iidx]) > 0){
                neighbors_x.push_back(neighborsX[idx][iidx]);
                if (neighbors_x.size() > max_r) break;
              }
            }
            if (neighbors_x.size() > n_excluded) {
              neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
            }

            // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_y;
            for (size_t iidx = 0; iidx < neighborsY[idx].size(); ++iidx){
              if(lib_set.count(neighborsY[idx][iidx])){
                neighbors_y.push_back(neighborsY[idx][iidx]);
                if (neighbors_y.size() > max_r) break;
              }
            }
            if (neighbors_y.size() > n_excluded) {
              neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
            }

            // Precompute y-neighbors set for fast lookup
            std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

            // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
            std::unordered_map<size_t, std::vector<size_t>> mapped_neighbors;

            for (size_t nx : neighbors_x) {
              if (neighborsY[nx][0] != std::numeric_limits<size_t>::max()) {
                for (size_t iidx = 0; iidx < neighborsY[nx].size(); ++iidx) {
                  if (lib_set.count(neighborsY[nx][iidx])) {
                    mapped_neighbors[nx].push_back(neighborsY[nx][iidx]);
                    if (mapped_neighbors[nx].size() >= num_neighbors) break;
                  }
                }
              }
            }

            // Compute intersection ratio between mapped x-neighbors and original y-neighbors
            for (size_t ki = 0; ki < num_neighbors; ++ki) {
              size_t count = 0;

              for (size_t nx : neighbors_x) {
                auto it = mapped_neighbors.find(nx);
                if (it != mapped_neighbors.end() && ki < it->second.size()) {
                  const auto& yn = it->second;

                  // Check if any of the first (ki+1) mapped neighbors exist in y's original neighbors
                  for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
                    if (y_neighbors_set.count(yn[pos])) {
                      ++count;
                      break;  // Count each x-neighbor only once if it intersects
                    }
                  }
                }
              }

              if (neighbors_x.size() > 0) {
                ratio_curves[i][ki] = static_cast<double>(count) / static_cast<double>(neighbors_x.size());
              }
            }
          }
        }
      }
    }

    std::vector<double> H1sequence;
    for (size_t col = 0; col < num_neighbors; ++col) {
      std::vector<double> mean_intersect;
      for (size_t row = 0; row < ratio_curves.size(); ++row){
        mean_intersect.push_back(ratio_curves[row][col]);
      }
      H1sequence.push_back(CppMean(mean_intersect,true));
    }

    return H1sequence;
  };

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<IntersectionRes> x_xmap_y;
    x_xmap_y.emplace_back(lib_size, ICSingle(lib_indices));
    return x_xmap_y;
  } else if (parallel_level == 1){
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<IntersectionRes> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      IntersectionRes result(lib_size, ICSingle(valid_lib_indices[i]));
      x_xmap_y[i] = result;
    }, threads);
    return x_xmap_y;
  } else {
    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<IntersectionRes> x_xmap_y;

    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Setup changing library
      if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }

      // Run cross map and store results
      x_xmap_y.emplace_back(lib_size, ICSingle(local_lib_indices));
    }

    return x_xmap_y;
  }
}

/*
 * Computes the Intersection Cardinality (IC) scores
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   lib: Library index vector (1-based in R, converted to 0-based).
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Number of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Number of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   - A vector representing the intersection cardinality (IC) scores, normalized between [0, 1].
 */
std::vector<double> IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int num_neighbors,
    int n_excluded,
    int threads,
    bool progressbar) {

  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return {0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  // Filter valid prediction points (exclude those with all NaN values)
  std::vector<int> valid_pred;
  for (int idx : pred) {
    if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;

    bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                             [](double v) { return std::isnan(v); });
    bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                             [](double v) { return std::isnan(v); });
    if (!x_nan && !y_nan) valid_pred.push_back(idx);
  }
  if (valid_pred.empty()) return {0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  // Parameter initialization
  const size_t k = static_cast<size_t>(num_neighbors);
  const size_t n_excluded_sizet = static_cast<size_t>(n_excluded);
  const size_t max_r = static_cast<size_t>(num_neighbors + n_excluded); // Total number of neighbors = actual used + excluded ones

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
  auto dist_x = CppMatDistance(embedding_x, false, true);
  auto dist_y = CppMatDistance(embedding_y, false, true);

  // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
  std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, std::numeric_limits<double>::quiet_NaN()));

  // Main parallel computation logic
  auto CMCSingle = [&](size_t i) {
    const int idx = valid_pred[i];

    // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
    auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r, lib);
    if (neighbors_x.size() > n_excluded_sizet) {
      neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
    }
    neighbors_x.resize(k); // Keep only the k actual neighbors

    // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
    auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r, lib);
    if (neighbors_y.size() > n_excluded_sizet) {
      neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
    }
    neighbors_y.resize(k); // Keep only the k actual neighbors

    // Precompute y-neighbors set for fast lookup
    std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

    // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
    std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
    for (size_t nx : neighbors_x) {
      mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k, lib);
    }

    // Compute intersection ratio between mapped x-neighbors and original y-neighbors
    for (size_t ki = 0; ki < k; ++ki) {
      size_t count = 0;
      for (size_t nx : neighbors_x) {
        if (ki < mapped_neighbors[nx].size()) {
          auto& yn = mapped_neighbors[nx];
          // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
          for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
            if (y_neighbors_set.count(yn[pos])) {
              ++count;
              break; // Count each x-neighbor only once if has any intersection
            }
          }
        }
      }
      if (!neighbors_x.empty()) {
        ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
      }
    }
  };

  if (progressbar) {
    // Parallel computation with a progress bar
    RcppThread::ProgressBar bar(valid_pred.size(), 1);
    RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
      CMCSingle(i);
      bar++;
    }, threads_sizet);
  } else {
    // Parallel computation without a progress bar
    RcppThread::parallelFor(0, valid_pred.size(), CMCSingle, threads_sizet);
  }

  std::vector<double> H1sequence;
  for (size_t col = 0; col < k; ++col) {
    std::vector<double> mean_intersect;
    for (size_t row = 0; row < ratio_curves.size(); ++row){
      mean_intersect.push_back(ratio_curves[row][col]);
    }
    H1sequence.push_back(CppMean(mean_intersect,true));
  }
  return H1sequence;
}

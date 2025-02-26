#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <unordered_set>
#include "CppStats.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

/*
 * Computes the Intersection Cardinality (IC) causal strength score.
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Number of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Number of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   - A vector containing the results of the DeLong test for the AUC values: [IC score, p-value, confidence interval upper bound, confidence interval lower bound].
 *   - IC score represents the causal strength score, normalized between [0, 1].
 *   - p-value indicates the statistical significance of the result.
 *   - The confidence interval bounds represent the uncertainty of the IC score.
 */
std::vector<double> IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
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
  size_t threads_sizet = static_cast<size_t>(threads);
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
  auto dist_x = CppMatDistance(embedding_x, false, true);
  auto dist_y = CppMatDistance(embedding_y, false, true);

  // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
  std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));

  // Main parallel computation logic
  auto CMCSingle = [&](size_t i) {
    const int idx = valid_pred[i];

    // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
    auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
    if (neighbors_x.size() > n_excluded_sizet) {
      neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
    }
    neighbors_x.resize(k); // Keep only the k actual neighbors

    // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
    auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
    if (neighbors_y.size() > n_excluded_sizet) {
      neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
    }
    neighbors_y.resize(k); // Keep only the k actual neighbors

    // Precompute y-neighbors set for fast lookup
    std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

    // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
    std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
    for (size_t nx : neighbors_x) {
      mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
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
      ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
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

  // Compute AUC (corresponding to ratio_to_auc in python package crossmapy)
  std::vector<double> H0sequence;
  for (size_t i = 1; i <= k; ++i) {
    H0sequence.push_back(static_cast<double>(i) / k);
  }

  std::vector<double> H1sequence;
  for (size_t col = 0; col < k; ++col) {
    std::vector<double> mean_intersect;
    for (size_t row = 0; row < ratio_curves.size(); ++row){
      mean_intersect.push_back(ratio_curves[row][col]);
    }
    H1sequence.push_back(CppMean(mean_intersect,true));
  }

  std::vector<double> dp_res = CppDeLongTest(H1sequence,H0sequence,">");

  return dp_res;
}

/**
 * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   A vector the results of the DeLong test for the AUC values: [IC score, p-value, confidence interval upper bound, confidence interval lower bound] one for each entry in num_neighbors.
 *   The result contains multiple rows, each corresponding to a different number of neighbors.
 */
std::vector<std::vector<double>> CrossMappingCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    const std::vector<int>& num_neighbors,
    const std::vector<int>& n_excluded,
    int threads,
    bool progressbar) {
  // Store results for each num_neighbors
  std::vector<std::vector<double>> results(num_neighbors.size(), std::vector<double>(4,std::numeric_limits<double>::quiet_NaN()));

  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return results;
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
  if (valid_pred.empty()) return results;

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(threads);
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
  auto dist_x = CppMatDistance(embedding_x, false, true);
  auto dist_y = CppMatDistance(embedding_y, false, true);

  auto CMCSingle = [&](size_t j) {
    const size_t k = static_cast<size_t>(num_neighbors[j]);
    const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
    const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones

    // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
    std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
      const int idx = valid_pred[i];

      // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
      auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
      if (neighbors_x.size() > n_excluded_sizet) {
        neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
      }
      neighbors_x.resize(k); // Keep only the k actual neighbors

      // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
      auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
      if (neighbors_y.size() > n_excluded_sizet) {
        neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
      }
      neighbors_y.resize(k); // Keep only the k actual neighbors

      // Precompute y-neighbors set for fast lookup
      std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

      // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
      std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
      for (size_t nx : neighbors_x) {
        mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
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
        ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
      }
    }, threads_sizet);

    std::vector<double> H0sequence;
    // for (size_t i = 1; i <= k; ++i) {
    //   H0sequence.push_back(static_cast<double>(i) / k);
    // }
    for (size_t i = 0; i < k; ++i) {
      H0sequence.push_back(static_cast<double>(i) / k);
    }

    std::vector<double> H1sequence;
    for (size_t col = 0; col < k; ++col) {
      std::vector<double> mean_intersect;
      for (size_t row = 0; row < ratio_curves.size(); ++row){
        mean_intersect.push_back(ratio_curves[row][col]);
      }
      H1sequence.push_back(CppMean(mean_intersect,true));
    }

    std::vector<double> dp_res = CppDeLongTest(H1sequence,H0sequence,">");
    results[j] = dp_res;
  };

  // Parallel computation with or without a progress bar
  if (progressbar) {
    RcppThread::ProgressBar bar(num_neighbors.size(), 1);
    for (size_t i = 0; i < num_neighbors.size(); ++i) {
      CMCSingle(i);
      bar++;
    };
  } else {
    for (size_t i = 0; i < num_neighbors.size(); ++i) {
      CMCSingle(i);
    }
  }

  return results; // Return the vector of results
}

/**
 * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
 *
 * Parameters:
 *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
 *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
 *   pred: Prediction index vector (1-based in R, converted to 0-based).
 *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
 *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
 *   threads: Number of parallel threads.
 *   progressbar: Whether to display a progress bar.
 *
 * Returns:
 *   A vector the results of the DeLong test for the AUC values: [IC score, p-value, confidence interval upper bound, confidence interval lower bound] one for each entry in num_neighbors.
 *   The result contains multiple rows, each corresponding to a different number of neighbors.
 */
std::vector<std::vector<double>> CrossMappingCardinality2(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    const std::vector<int>& num_neighbors,
    const std::vector<int>& n_excluded,
    int threads,
    bool progressbar) {
  // Store results for each num_neighbors
  std::vector<std::vector<double>> results(num_neighbors.size(), std::vector<double>(4,std::numeric_limits<double>::quiet_NaN()));

  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return results;
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
  if (valid_pred.empty()) return results;

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(threads);
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
  auto dist_x = CppMatDistance(embedding_x, false, true);
  auto dist_y = CppMatDistance(embedding_y, false, true);

  // Main parallel computation logic for each (num_neighbors, n_excluded) pair
  auto CMCSingle = [&](size_t j) {
    const size_t k = static_cast<size_t>(num_neighbors[j]);
    const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
    const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones

    // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
    std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));

    for (size_t i = 0; i < valid_pred.size(); ++i) {
      const int idx = valid_pred[i];

      // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
      auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
      if (neighbors_x.size() > n_excluded_sizet) {
        neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
      }
      neighbors_x.resize(k); // Keep only the k actual neighbors

      // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
      auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
      if (neighbors_y.size() > n_excluded_sizet) {
        neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
      }
      neighbors_y.resize(k); // Keep only the k actual neighbors

      // Precompute y-neighbors set for fast lookup
      std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

      // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
      std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
      for (size_t nx : neighbors_x) {
        mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
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
        ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
      }
    }

    // Compute AUC for the current (num_neighbors, n_excluded) pair
    std::vector<double> H0sequence;
    // for (size_t i = 1; i <= k; ++i) {
    //   H0sequence.push_back(static_cast<double>(i) / k);
    // }
    for (size_t i = 0; i < k; ++i) {
      H0sequence.push_back(static_cast<double>(i) / k);
    }

    std::vector<double> H1sequence;
    for (size_t col = 0; col < k; ++col) {
      std::vector<double> mean_intersect;
      for (size_t row = 0; row < ratio_curves.size(); ++row){
        mean_intersect.push_back(ratio_curves[row][col]);
      }
      H1sequence.push_back(CppMean(mean_intersect,true));
    }

    std::vector<double> dp_res = CppDeLongTest(H1sequence,H0sequence,">");
    results[j] = dp_res;
  };

  // Parallel computation with or without a progress bar
  if (progressbar) {
    RcppThread::ProgressBar bar(num_neighbors.size(), 1);
    RcppThread::parallelFor(0, num_neighbors.size(), [&](size_t j) {
      CMCSingle(j);
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, num_neighbors.size(), CMCSingle, threads_sizet);
  }

  return results; // Return the vector of results
}

// The previous CrossMappingCardinality implementation is retained for comparison.
//    ----- Wenbo Lv, 2025.02.25

// #include <vector>
// #include <cmath>
// #include <algorithm>
// #include <numeric>
// #include <limits>
// #include <utility>
// #include <unordered_set>
// #include "CppStats.h"
// #include <RcppThread.h>
//
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(RcppThread)]]
//
// /*
//  * Computes the Intersection Cardinality (IC) causal strength score.
//  *
//  * Parameters:
//  *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
//  *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
//  *   pred: Prediction index vector (1-based in R, converted to 0-based).
//  *   num_neighbors: Number of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
//  *   n_excluded: Number of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
//  *   threads: Number of parallel threads.
//  *   progressbar: Whether to display a progress bar.
//  *
//  * Returns:
//  *   - A double representing the IC causal strength score, normalized between [0,1].
//  */
// double IntersectionCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& pred,
//     int num_neighbors,
//     int n_excluded,
//     int threads,
//     bool progressbar) {
//
//   // Input validation
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return 0.0;
//   }
//
//   // Filter valid prediction points (exclude those with all NaN values)
//   std::vector<int> valid_pred;
//   for (int idx : pred) {
//     if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;
//
//     bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     if (!x_nan && !y_nan) valid_pred.push_back(idx);
//   }
//   if (valid_pred.empty()) return 0.0;
//
//   // Parameter initialization
//   const size_t k = static_cast<size_t>(num_neighbors);
//   const size_t n_excluded_sizet = static_cast<size_t>(n_excluded);
//   const size_t max_r = static_cast<size_t>(num_neighbors + n_excluded); // Total number of neighbors = actual used + excluded ones
//
//   // Configure threads
//   size_t threads_sizet = static_cast<size_t>(threads);
//   threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);
//
//   // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
//   auto dist_x = CppMatDistance(embedding_x, false, true);
//   auto dist_y = CppMatDistance(embedding_y, false, true);
//
//   // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//   std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));
//
//   // Main parallel computation logic
//   auto CMCSingle = [&](size_t i) {
//     const int idx = valid_pred[i];
//
//     // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//     auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
//     if (neighbors_x.size() > n_excluded_sizet) {
//       neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
//     }
//     neighbors_x.resize(k); // Keep only the k actual neighbors
//
//     // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//     auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
//     if (neighbors_y.size() > n_excluded_sizet) {
//       neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
//     }
//     neighbors_y.resize(k); // Keep only the k actual neighbors
//
//     // Precompute y-neighbors set for fast lookup
//     std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//
//     // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//     std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//     for (size_t nx : neighbors_x) {
//       mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
//     }
//
//     // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//     for (size_t ki = 0; ki < k; ++ki) {
//       size_t count = 0;
//       for (size_t nx : neighbors_x) {
//         if (ki < mapped_neighbors[nx].size()) {
//           auto& yn = mapped_neighbors[nx];
//           // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//           for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//             if (y_neighbors_set.count(yn[pos])) {
//               ++count;
//               break; // Count each x-neighbor only once if has any intersection
//             }
//           }
//         }
//       }
//       ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//     }
//   };
//
//   // auto CMCSingle = [&](size_t i) {
//   //   const int idx = valid_pred[i];
//   //
//   //   // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//   //   auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
//   //   if (neighbors_x.size() > n_excluded_sizet) {
//   //     neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
//   //   }
//   //   neighbors_x.resize(k); // Keep only the k actual neighbors
//   //
//   //   // Retrieve y's neighbor indices (for mapping validation)
//   //   std::vector<std::vector<size_t>> mapped_neighbors(embedding_y.size());
//   //   for (size_t nx : neighbors_x) {
//   //     mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
//   //   }
//   //
//   //   // Compute mapping ratio for each k (corresponding to count_mapping in python package crossmapy)
//   //   for (size_t ki = 0; ki < k; ++ki) {
//   //     size_t count = 0;
//   //     for (size_t nx : neighbors_x) {
//   //       if (ki < mapped_neighbors[nx].size()) {
//   //         auto& yn = mapped_neighbors[nx];
//   //         if (std::find(yn.begin(), yn.begin() + ki + 1, idx) != yn.begin() + ki + 1) {
//   //           ++count;
//   //         }
//   //       }
//   //     }
//   //     ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//   //   }
//   // };
//
//   if (progressbar) {
//     // Parallel computation with a progress bar
//     RcppThread::ProgressBar bar(valid_pred.size(), 1);
//     RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
//       CMCSingle(i);
//       bar++;
//     }, threads_sizet);
//   } else {
//     // Parallel computation without a progress bar
//     RcppThread::parallelFor(0, valid_pred.size(), CMCSingle, threads_sizet);
//   }
//
//   // Compute AUC (corresponding to ratio_to_auc in python package crossmapy)
//   double total_auc = 0.0;
//   for (const auto& curve : ratio_curves) {
//     double auc = 0.0;
//     for (size_t i = 1; i < curve.size(); ++i) {
//       auc += (curve[i] + curve[i - 1]) * 0.5;
//     }
//     auc /= (curve.size() - 1); // Normalize
//     total_auc += auc;
//   }
//   double mean_auc = total_auc / ratio_curves.size();
//
//   // Convert to final score (corresponding to auc_to_score in python package crossmapy)
//   // double cmc_score = std::max(0.0, mean_auc); // just to ensure that the causal score is a number greater than 0.
//   double cmc_score = 2.0 * (mean_auc - 0.5);
//   // if (cmc_score < 0) {
//   //   cmc_score = 2.0 * (0.5 - mean_auc);
//   // }
//
//   return std::max(0.0, cmc_score); // Ensure non-negative result
// }
//
// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
//  *
//  * Parameters:
//  *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
//  *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
//  *   pred: Prediction index vector (1-based in R, converted to 0-based).
//  *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
//  *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
//  *   threads: Number of parallel threads.
//  *   progressbar: Whether to display a progress bar.
//  *
//  * Returns:
//  *   A vector of normalized CMC causal strength scores in the range [0,1], corresponding to each value in num_neighbors.
//  */
// std::vector<double> CrossMappingCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& pred,
//     const std::vector<int>& num_neighbors,
//     const std::vector<int>& n_excluded,
//     int threads,
//     bool progressbar) {
//   // Store results for each num_neighbors
//   std::vector<double> results(num_neighbors.size(), 0.0);
//
//   // Input validation
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return results;
//   }
//
//   // Filter valid prediction points (exclude those with all NaN values)
//   std::vector<int> valid_pred;
//   for (int idx : pred) {
//     if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;
//
//     bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     if (!x_nan && !y_nan) valid_pred.push_back(idx);
//   }
//   if (valid_pred.empty()) return results;
//
//   // Configure threads
//   size_t threads_sizet = static_cast<size_t>(threads);
//   threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);
//
//   // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
//   auto dist_x = CppMatDistance(embedding_x, false, true);
//   auto dist_y = CppMatDistance(embedding_y, false, true);
//
//   // Main parallel computation logic for each (num_neighbors, n_excluded) pair
//   // auto CMCSingle = [&](size_t j) {
//   //   const size_t k = static_cast<size_t>(num_neighbors[j]);
//   //   const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
//   //   const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones
//   //
//   //   // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//   //   std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));
//   //
//   //   for (size_t i = 0; i < valid_pred.size(); ++i) {
//   //     const int idx = valid_pred[i];
//   //
//   //     // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//   //     auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
//   //     if (neighbors_x.size() > n_excluded_sizet) {
//   //       neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
//   //     }
//   //     neighbors_x.resize(k); // Keep only the k actual neighbors
//   //
//   //     // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//   //     auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
//   //     if (neighbors_y.size() > n_excluded_sizet) {
//   //       neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
//   //     }
//   //     neighbors_y.resize(k); // Keep only the k actual neighbors
//   //
//   //     // Precompute y-neighbors set for fast lookup
//   //     std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//   //
//   //     // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//   //     std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//   //     for (size_t nx : neighbors_x) {
//   //       mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
//   //     }
//   //
//   //     // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//   //     for (size_t ki = 0; ki < k; ++ki) {
//   //       size_t count = 0;
//   //       for (size_t nx : neighbors_x) {
//   //         if (ki < mapped_neighbors[nx].size()) {
//   //           auto& yn = mapped_neighbors[nx];
//   //           // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//   //           for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//   //             if (y_neighbors_set.count(yn[pos])) {
//   //               ++count;
//   //               break; // Count each x-neighbor only once if has any intersection
//   //             }
//   //           }
//   //         }
//   //       }
//   //       ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//   //     }
//   //   }
//   //
//   //   // Compute AUC for the current (num_neighbors, n_excluded) pair
//   //   double total_auc = 0.0;
//   //   for (const auto& curve : ratio_curves) {
//   //     double auc = 0.0;
//   //     for (size_t i = 1; i < curve.size(); ++i) {
//   //       auc += (curve[i] + curve[i - 1]) * 0.5;
//   //     }
//   //     auc /= (curve.size() - 1); // Normalize
//   //     total_auc += auc;
//   //   }
//   //   double mean_auc = total_auc / ratio_curves.size();
//   //
//   //   // Convert to final score (corresponding to auc_to_score in python package crossmapy)
//   //   double cmc_score = 2.0 * (mean_auc - 0.5);
//   //   results[j] = std::max(0.0, cmc_score); // Ensure non-negative result
//   // };
//
//   auto CMCSingle = [&](size_t j) {
//     const size_t k = static_cast<size_t>(num_neighbors[j]);
//     const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
//     const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones
//
//     // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//     std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));
//
//     // Perform the operations using RcppThread
//     RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
//       const int idx = valid_pred[i];
//
//       // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//       auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
//       if (neighbors_x.size() > n_excluded_sizet) {
//         neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
//       }
//       neighbors_x.resize(k); // Keep only the k actual neighbors
//
//       // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//       auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
//       if (neighbors_y.size() > n_excluded_sizet) {
//         neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
//       }
//       neighbors_y.resize(k); // Keep only the k actual neighbors
//
//       // Precompute y-neighbors set for fast lookup
//       std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//
//       // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//       std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//       for (size_t nx : neighbors_x) {
//         mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
//       }
//
//       // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//       for (size_t ki = 0; ki < k; ++ki) {
//         size_t count = 0;
//         for (size_t nx : neighbors_x) {
//           if (ki < mapped_neighbors[nx].size()) {
//             auto& yn = mapped_neighbors[nx];
//             // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//             for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//               if (y_neighbors_set.count(yn[pos])) {
//                 ++count;
//                 break; // Count each x-neighbor only once if has any intersection
//               }
//             }
//           }
//         }
//         ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//       }
//     }, threads_sizet);
//
//     // Compute AUC for the current (num_neighbors, n_excluded) pair
//     double total_auc = 0.0;
//     for (const auto& curve : ratio_curves) {
//       double auc = 0.0;
//       for (size_t i = 1; i < curve.size(); ++i) {
//         auc += (curve[i] + curve[i - 1]) * 0.5;
//       }
//       auc /= (curve.size() - 1); // Normalize
//       total_auc += auc;
//     }
//     double mean_auc = total_auc / ratio_curves.size();
//
//     // Convert to final score (corresponding to auc_to_score in python package crossmapy)
//     double cmc_score = 2.0 * (mean_auc - 0.5);
//     results[j] = std::max(0.0, cmc_score); // Ensure non-negative result
//   };
//
//   // Parallel computation with or without a progress bar
//   if (progressbar) {
//     RcppThread::ProgressBar bar(num_neighbors.size(), 1);
//     for (size_t i = 0; i < num_neighbors.size(); ++i) {
//       CMCSingle(i);
//       bar++;
//     };
//   } else {
//     for (size_t i = 0; i < num_neighbors.size(); ++i) {
//       CMCSingle(i);
//     }
//   }
//
//   return results; // Return the vector of results
// }
//
// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
//  *
//  * Parameters:
//  *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
//  *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
//  *   pred: Prediction index vector (1-based in R, converted to 0-based).
//  *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
//  *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
//  *   threads: Number of parallel threads.
//  *   progressbar: Whether to display a progress bar.
//  *
//  * Returns:
//  *   A vector of normalized CMC causal strength scores in the range [0,1], corresponding to each value in num_neighbors.
//  */
// std::vector<double> CrossMappingCardinality2(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& pred,
//     const std::vector<int>& num_neighbors,
//     const std::vector<int>& n_excluded,
//     int threads,
//     bool progressbar) {
//   // Store results for each num_neighbors
//   std::vector<double> results(num_neighbors.size(), 0.0);
//
//   // Input validation
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return results;
//   }
//
//   // Filter valid prediction points (exclude those with all NaN values)
//   std::vector<int> valid_pred;
//   for (int idx : pred) {
//     if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;
//
//     bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     if (!x_nan && !y_nan) valid_pred.push_back(idx);
//   }
//   if (valid_pred.empty()) return results;
//
//   // Configure threads
//   size_t threads_sizet = static_cast<size_t>(threads);
//   threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);
//
//   // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
//   auto dist_x = CppMatDistance(embedding_x, false, true);
//   auto dist_y = CppMatDistance(embedding_y, false, true);
//
//   // Main parallel computation logic for each (num_neighbors, n_excluded) pair
//   auto CMCSingle = [&](size_t j) {
//     const size_t k = static_cast<size_t>(num_neighbors[j]);
//     const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
//     const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones
//
//     // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//     std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, 0.0));
//
//     for (size_t i = 0; i < valid_pred.size(); ++i) {
//       const int idx = valid_pred[i];
//
//       // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//       auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
//       if (neighbors_x.size() > n_excluded_sizet) {
//         neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
//       }
//       neighbors_x.resize(k); // Keep only the k actual neighbors
//
//       // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//       auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r);
//       if (neighbors_y.size() > n_excluded_sizet) {
//         neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
//       }
//       neighbors_y.resize(k); // Keep only the k actual neighbors
//
//       // Precompute y-neighbors set for fast lookup
//       std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//
//       // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//       std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//       for (size_t nx : neighbors_x) {
//         mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
//       }
//
//       // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//       for (size_t ki = 0; ki < k; ++ki) {
//         size_t count = 0;
//         for (size_t nx : neighbors_x) {
//           if (ki < mapped_neighbors[nx].size()) {
//             auto& yn = mapped_neighbors[nx];
//             // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//             for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//               if (y_neighbors_set.count(yn[pos])) {
//                 ++count;
//                 break; // Count each x-neighbor only once if has any intersection
//               }
//             }
//           }
//         }
//         ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//       }
//     }
//
//     // Compute AUC for the current (num_neighbors, n_excluded) pair
//     double total_auc = 0.0;
//     for (const auto& curve : ratio_curves) {
//       double auc = 0.0;
//       for (size_t i = 1; i < curve.size(); ++i) {
//         auc += (curve[i] + curve[i - 1]) * 0.5;
//       }
//       auc /= (curve.size() - 1); // Normalize
//       total_auc += auc;
//     }
//     double mean_auc = total_auc / ratio_curves.size();
//
//     // Convert to final score (corresponding to auc_to_score in python package crossmapy)
//     double cmc_score = 2.0 * (mean_auc - 0.5);
//     results[j] = std::max(0.0, cmc_score); // Ensure non-negative result
//   };
//
//   // Parallel computation with or without a progress bar
//   if (progressbar) {
//     RcppThread::ProgressBar bar(num_neighbors.size(), 1);
//     RcppThread::parallelFor(0, num_neighbors.size(), [&](size_t j) {
//       CMCSingle(j);
//       bar++;
//     }, threads_sizet);
//   } else {
//     RcppThread::parallelFor(0, num_neighbors.size(), CMCSingle, threads_sizet);
//   }
//
//   return results; // Return the vector of results
// }

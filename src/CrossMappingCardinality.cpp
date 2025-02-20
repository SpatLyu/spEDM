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

/**
 * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
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
 *   Normalized CMC causal strength score in the range [0,1].
 */
double CrossMappingCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    int num_neighbors,
    int n_excluded,
    int threads,
    bool progressbar) {

  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return 0.0;
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
  if (valid_pred.empty()) return 0.0;

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

  // auto CMCSingle = [&](size_t i) {
  //   const int idx = valid_pred[i];
  //
  //   // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
  //   auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r);
  //   if (neighbors_x.size() > n_excluded_sizet) {
  //     neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
  //   }
  //   neighbors_x.resize(k); // Keep only the k actual neighbors
  //
  //   // Retrieve y's neighbor indices (for mapping validation)
  //   std::vector<std::vector<size_t>> y_neighbors(embedding_y.size());
  //   for (size_t nx : neighbors_x) {
  //     y_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k);
  //   }
  //
  //   // Compute mapping ratio for each k (corresponding to count_mapping in python package crossmapy)
  //   for (size_t ki = 0; ki < k; ++ki) {
  //     size_t count = 0;
  //     for (size_t nx : neighbors_x) {
  //       if (ki < y_neighbors[nx].size()) {
  //         auto& yn = y_neighbors[nx];
  //         if (std::find(yn.begin(), yn.begin() + ki + 1, idx) != yn.begin() + ki + 1) {
  //           ++count;
  //         }
  //       }
  //     }
  //     ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
  //   }
  // };

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
  double total_auc = 0.0;
  for (const auto& curve : ratio_curves) {
    double auc = 0.0;
    for (size_t i = 1; i < curve.size(); ++i) {
      auc += (curve[i] + curve[i - 1]) * 0.5;
    }
    auc /= (curve.size() - 1); // Normalize
    total_auc += auc;
  }
  double mean_auc = total_auc / ratio_curves.size();

  // Convert to final score (corresponding to auc_to_score in python package crossmapy)
  double cmc_score = 2.0 * (mean_auc - 0.5);
  // double cmc_score = std::max(0.0, mean_auc);
  // if (cmc_score < 0) {
  //   cmc_score = 2.0 * (0.5 - mean_auc);
  // }
  return std::max(0.0, cmc_score); // Ensure non-negative result
}

/*
 * Computes the Intersection Cardinality (IC) causal strength score.
 *
 * Parameters:
 *   - embedding_x:   The state-space reconstructed from the potential cause variable.
 *   - embedding_y:   The state-space reconstructed from the potential effect variable.
 *   - pred:          A vector specifying the prediction indices(1-based in R, converted to 0-based in C++).
 *   - num_neighbors: Number of neighbors used for cross-mapping.
 *   - max_neighbors: Maximum number of neighbors usable for IC computation.
 *   - threads:       Number of threads to use for parallel processing.
 *   - progressbar:   If true, display a progress bar during computation.
 *
 * Returns:
 *   - A double representing the IC causal strength score, normalized between [0,1].
 */
double IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<int>& pred,
    int num_neighbors,
    int max_neighbors,
    int threads,
    bool progressbar
) {
  // Construct a valid_pred vector to store indices that are not entirely NaN
  std::vector<int> valid_pred;

  for (int idx : pred) {
    // Ensure index is within valid range
    if (idx < 0 || static_cast<std::size_t>(idx) >= embedding_x.size()) {
      continue;
    }

    // Check if all values in embedding_x[idx] and embedding_y[idx] are NaN
    bool x_all_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
    bool y_all_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });

    // If at least one of them has valid data, add idx to valid_pred
    if (!x_all_nan && !y_all_nan) {
      valid_pred.push_back(idx);
    }
  }

  // If no valid predictions remain, return 0.0
  if (valid_pred.empty()) {
    return 0.0;
  }

  std::size_t k = static_cast<size_t>(num_neighbors);
  std::size_t max_r = std::min(static_cast<size_t>(max_neighbors), embedding_x.size());
  max_r = std::max(max_r,k);

  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  // Compute the Intersectional Cardinality (IC) curve
  std::vector<std::pair<int, double>> IC_curve;

  // // Iterate over each r
  // for (size_t r = k; r <= max_r; ++r) {
  //   double intersection_sum = 0.0;
  //
  //   for (size_t t = 0; t < pred.size(); ++t) {
  //     // Find k-nearest neighbors in embedding_x and embedding_y
  //     std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, pred[t]-1, k);
  //     std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, pred[t]-1, r);
  //
  //     // For each neighbor in embedding_x, find its corresponding neighbor in embedding_y
  //     std::vector<std::size_t> neighbors_xmapped;
  //     for (std::size_t nx : neighbors_x) {
  //       std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, nx, 1); // Map to 1 nearest neighbor
  //       if (!neighbors_y.empty()) {
  //         neighbors_xmapped.push_back(neighbors_y[0]);
  //       }
  //     }
  //
  //     // Compute the intersection count
  //     double intersection_count = 0;
  //     for (size_t nxm : neighbors_xmapped) {
  //       if (std::find(neighbors_y.begin(), neighbors_y.end(), nxm) != neighbors_y.end()) {
  //         intersection_count++;
  //       }
  //     }
  //
  //     intersection_sum += intersection_count;
  //   }
  //
  //   IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
  // }

  // Perform the operations using RcppThread
  if (progressbar) {
    RcppThread::ProgressBar bar(max_neighbors - num_neighbors - 1, 1);
    RcppThread::parallelFor(k, max_r, [&](size_t r) {
      double intersection_sum = 0.0;

      for (size_t t = 0; t < pred.size(); ++t) {
        // Find k-nearest neighbors in embedding_x and embedding_y
        std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, pred[t]-1, k);
        std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, pred[t]-1, r);

        // For each neighbor in embedding_x, find its corresponding neighbor in embedding_y
        std::vector<std::size_t> neighbors_xmapped;
        for (std::size_t nx : neighbors_x) {
          std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, nx, 1); // Map to 1 nearest neighbor
          if (!neighbors_y.empty()) {
            neighbors_xmapped.push_back(neighbors_y[0]);
          }
        }

        // Compute the intersection count
        double intersection_count = 0;
        for (size_t nxm : neighbors_xmapped) {
          if (std::find(neighbors_y.begin(), neighbors_y.end(), nxm) != neighbors_y.end()) {
            intersection_count++;
          }
        }

        intersection_sum += intersection_count;
      }

      IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(k, max_r, [&](size_t r) {
      double intersection_sum = 0.0;

      for (size_t t = 0; t < pred.size(); ++t) {
        // Find k-nearest neighbors in embedding_x and embedding_y
        std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, pred[t]-1, k);
        std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, pred[t]-1, r);

        // For each neighbor in embedding_x, find its corresponding neighbor in embedding_y
        std::vector<std::size_t> neighbors_xmapped;
        for (std::size_t nx : neighbors_x) {
          std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, nx, 1); // Map to 1 nearest neighbor
          if (!neighbors_y.empty()) {
            neighbors_xmapped.push_back(neighbors_y[0]);
          }
        }

        // Compute the intersection count
        double intersection_count = 0;
        for (size_t nxm : neighbors_xmapped) {
          if (std::find(neighbors_y.begin(), neighbors_y.end(), nxm) != neighbors_y.end()) {
            intersection_count++;
          }
        }

        intersection_sum += intersection_count;
      }

      IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
    }, threads_sizet);
  }

  // Sort the vector based on the first element (int) in ascending order
  std::sort(IC_curve.begin(), IC_curve.end(),
            [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
              return a.first < b.first;
            });

  // Remove all pairs where the double value is >= 1
  IC_curve.erase(std::remove_if(IC_curve.begin(), IC_curve.end(),
                                [](const std::pair<int, double>& p) {
                                  return p.second >= 1;  // Remove if double is >= 1
                                }), IC_curve.end());

  // Take out IC value separately
  std::vector<double> IC_Values;
  for (const auto& result : IC_curve) {
    IC_Values.push_back(result.second);
  }

  // Counter to track occurrences of 1.0
  int count = 0;

  // Remove only the second occurrence of 1.0
  auto it = std::remove_if(IC_Values.begin(), IC_Values.end(),
                           [&count](double value) {
                             if (value == 1.0) {
                               count++;
                               return count == 2;  // Remove only the second occurrence
                             }
                             return false;
                           });

  // Erase the marked element
  if (it != IC_Values.end()) {
    IC_Values.erase(it, IC_Values.end());
  }

  // Compute AUC and normalize causal strength score
  double auc = 0.0;
  double r_size = static_cast<double>(IC_Values.size()) - 1;

  for (size_t i = 1; i < IC_Values.size(); ++i) {
    double width = 1.0;  // Assuming uniform step size in r
    auc += (IC_Values[i] + IC_Values[i - 1]) * width / 2.0;
  }

  double AUC = auc / r_size;
  double ICC_score = std::max(0.0, 2.0 * (AUC - 0.5));

  return ICC_score;
}

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
// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score.
//  *
//  * Parameters:
//  *   - embedding_x:   The state-space reconstructed from the potential cause variable.
//  *   - embedding_y:   The state-space reconstructed from the potential effect variable.
//  *   - pred:          A vector specifying the prediction indices(1-based in R, converted to 0-based in C++).
//  *   - num_neighbors: Number of neighbors used for cross-mapping.
//  *   - n_excluded:    Number of excluded neighbors in the distance matrix.
//  *   - threads:       Number of threads to use for parallel processing.
//  *   - progressbar:   If true, display a progress bar during computation.
//  *
//  * Returns:
//  *   - A double representing the CMC causal strength score, normalized between [0,1].
//  */
// double CrossMappingCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& pred,
//     int num_neighbors,
//     int n_excluded,
//     int threads,
//     bool progressbar
// ) {
//   // Ensure valid input dimensions
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return 0.0;
//   }
//
//   // Construct a valid_pred vector to store indices that are not entirely NaN
//   std::vector<int> valid_pred;
//
//   for (int idx : pred) {
//     // Ensure index is within valid range
//     if (idx < 0 || static_cast<std::size_t>(idx) >= embedding_x.size()) {
//       continue;
//     }
//
//     // Check if all values in embedding_x[idx] and embedding_y[idx] are NaN
//     bool x_all_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
//                                  [](double v) { return std::isnan(v); });
//     bool y_all_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
//                                  [](double v) { return std::isnan(v); });
//
//     if (!x_all_nan && !y_all_nan) {
//       valid_pred.push_back(idx);
//     }
//
//     // // If all of them have no nan value, add idx to valid_pred
//     // if (!checkOneDimVectorHasNaN(embedding_x[idx]) && !checkOneDimVectorHasNaN(embedding_y[idx])) {
//     //   valid_pred.push_back(idx);
//     // }
//   }
//
//   // If no valid predictions remain, return 0.0
//   if (valid_pred.empty()) {
//     return 0.0;
//   }
//
//   std::size_t k = static_cast<size_t>(num_neighbors);
//   // std::size_t max_r = std::min(static_cast<size_t>(num_neighbors+n_excluded), embedding_x.size());
//   std::size_t max_r = std::min(static_cast<size_t>(n_excluded), embedding_x.size());
//   max_r = std::max(max_r, static_cast<size_t>(num_neighbors + 1));
//
//   size_t threads_sizet = static_cast<size_t>(threads);
//   unsigned int max_threads = std::thread::hardware_concurrency();
//   threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);
//
//   // Compute distance matrices for embedding_x and embedding_y
//   std::vector<std::vector<double>> dist_x = CppMatDistance(embedding_x,false,true);
//   std::vector<std::vector<double>> dist_y = CppMatDistance(embedding_y,false,true);
//
//   // Compute the mapping ratios
//   std::vector<double> mapping_ratios(valid_pred.size(), 0.0);
//
//   // // Iterate over each valid_pred
//   // for (size_t i = 0; i < valid_pred.size(); ++i) {
//   //   std::vector<std::size_t> neighbors_x = CppDistKNNIndice(dist_x, valid_pred[i], k);
//   //   std::vector<std::size_t> neighbors_y = CppDistKNNIndice(dist_y, valid_pred[i], max_r);
//   //
//   //   // Map neighbors from embedding_x to embedding_y
//   //   std::unordered_set<size_t> mapped_neighbors;
//   //   for (size_t nx : neighbors_x) {
//   //     auto mapped = CppDistKNNIndice(dist_y, nx, 1);
//   //     if (!mapped.empty()) {
//   //       mapped_neighbors.insert(mapped[0]);
//   //     }
//   //   }
//   //
//   //   // Compute the intersection count
//   //   double intersection_count = 0;
//   //   for (size_t ny : neighbors_y) {
//   //     if (mapped_neighbors.find(ny) != mapped_neighbors.end()) {
//   //       intersection_count++;
//   //     }
//   //   }
//   //
//   //   mapping_ratios.emplace_back(intersection_count / neighbors_y.size());
//   // }
//
//   // Perform the operations using RcppThread
//   if (progressbar) {
//     RcppThread::ProgressBar bar(valid_pred.size(), 1);
//     RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
//       std::vector<std::size_t> neighbors_x = CppDistKNNIndice(dist_x, valid_pred[i], k);
//       std::vector<std::size_t> neighbors_y = CppDistKNNIndice(dist_y, valid_pred[i], max_r);
//
//       // Map neighbors from embedding_x to embedding_y
//       std::unordered_set<size_t> mapped_neighbors;
//       for (size_t nx : neighbors_x) {
//         auto mapped = CppDistKNNIndice(dist_y, nx, 1);
//         if (!mapped.empty()) {
//           mapped_neighbors.insert(mapped[0]);
//         }
//       }
//
//       // Compute the intersection count
//       double intersection_count = 0;
//       for (size_t ny : neighbors_y) {
//         if (mapped_neighbors.find(ny) != mapped_neighbors.end()) {
//           intersection_count++;
//         }
//       }
//
//       mapping_ratios[i] = intersection_count / neighbors_y.size();
//       bar++;
//     }, threads_sizet);
//   } else {
//     RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
//       std::vector<std::size_t> neighbors_x = CppDistKNNIndice(dist_x, valid_pred[i], k);
//       std::vector<std::size_t> neighbors_y = CppDistKNNIndice(dist_y, valid_pred[i], max_r);
//
//       // Map neighbors from embedding_x to embedding_y
//       std::unordered_set<size_t> mapped_neighbors;
//       for (size_t nx : neighbors_x) {
//         auto mapped = CppDistKNNIndice(dist_y, nx, 1);
//         if (!mapped.empty()) {
//           mapped_neighbors.insert(mapped[0]);
//         }
//       }
//
//       // Compute the intersection count
//       double intersection_count = 0;
//       for (size_t ny : neighbors_y) {
//         if (mapped_neighbors.find(ny) != mapped_neighbors.end()) {
//           intersection_count++;
//         }
//       }
//
//       mapping_ratios[i] = intersection_count / neighbors_y.size();
//     }, threads_sizet);
//   }
//
//   // Compute AUC and normalize the score
//   std::sort(mapping_ratios.begin(), mapping_ratios.end());
//   double auc = 0.0;
//   for (size_t i = 1; i < mapping_ratios.size(); ++i) {
//     auc += (mapping_ratios[i] + mapping_ratios[i - 1]) * 0.5;
//   }
//   auc /= mapping_ratios.size();
//
//   double cmc_score;
//   // Normalize the score to [0,1]
//   cmc_score = std::max(0.0, 2.0 * (auc - 0.5));
//   // cmc_score = std::max(0.0, auc);
//   if (cmc_score < 0) {
//    cmc_score = 2.0 * (0.5 - auc);
//   }
//
//   return cmc_score;
// }
//
// /*
//  * Computes the Intersection Cardinality (IC) causal strength score.
//  *
//  * Parameters:
//  *   - embedding_x:   The state-space reconstructed from the potential cause variable.
//  *   - embedding_y:   The state-space reconstructed from the potential effect variable.
//  *   - num_neighbors: Number of neighbors used for cross-mapping.
//  *   - max_neighbors: Maximum number of neighbors usable for IC computation.
//  *   - threads:       Number of threads to use for parallel processing.
//  *   - progressbar:   If true, display a progress bar during computation.
//  *
//  * Returns:
//  *   - A double representing the IC causal strength score, normalized between [0,1].
//  */
// double IntersectionCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     int num_neighbors,
//     int max_neighbors,
//     int threads,
//     bool progressbar
// ) {
//   std::size_t k = static_cast<size_t>(num_neighbors);
//   std::size_t max_r = static_cast<size_t>(max_neighbors);
//
//   size_t threads_sizet = static_cast<size_t>(threads);
//   unsigned int max_threads = std::thread::hardware_concurrency();
//   threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);
//
//   // Compute the Intersectional Cardinality (IC) curve
//   std::vector<std::pair<int, double>> IC_curve;
//
//   // // Iterate over each r
//   // for (size_t r = k; r <= max_r; ++r) {
//   //   double intersection_sum = 0.0;
//   //
//   //   for (size_t t = 0; t < embedding_x.size(); ++t) {
//   //     std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, t, k);
//   //     std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, t, r);
//   //
//   //     // Compute the intersection count
//   //     double intersection_count = 0;
//   //     for (size_t nx : neighbors_x) {
//   //       if (std::find(neighbors_y.begin(), neighbors_y.end(), nx) != neighbors_y.end()) {
//   //         intersection_count++;
//   //       }
//   //     }
//   //
//   //     intersection_sum += intersection_count;
//   //   }
//   //
//   //   IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
//   // }
//
//   // Perform the operations using RcppThread
//   if (progressbar) {
//     RcppThread::ProgressBar bar(max_neighbors - num_neighbors - 1, 1);
//     RcppThread::parallelFor(k, max_r, [&](size_t r) {
//       double intersection_sum = 0.0;
//
//       for (size_t t = 0; t < embedding_x.size(); ++t) {
//         std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, t, k);
//         std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, t, r);
//
//         // Compute the intersection count
//         double intersection_count = 0;
//         for (size_t nx : neighbors_x) {
//           if (std::find(neighbors_y.begin(), neighbors_y.end(), nx) != neighbors_y.end()) {
//             intersection_count++;
//           }
//         }
//
//         intersection_sum += intersection_count;
//       }
//
//       IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
//       bar++;
//     }, threads_sizet);
//   } else {
//     RcppThread::parallelFor(k, max_r, [&](size_t r) {
//       double intersection_sum = 0.0;
//
//       for (size_t t = 0; t < embedding_x.size(); ++t) {
//         std::vector<std::size_t> neighbors_x = CppKNNIndice(embedding_x, t, k);
//         std::vector<std::size_t> neighbors_y = CppKNNIndice(embedding_y, t, r);
//
//         // Compute the intersection count
//         double intersection_count = 0;
//         for (size_t nx : neighbors_x) {
//           if (std::find(neighbors_y.begin(), neighbors_y.end(), nx) != neighbors_y.end()) {
//             intersection_count++;
//           }
//         }
//
//         intersection_sum += intersection_count;
//       }
//
//       IC_curve.emplace_back(static_cast<int>(r), intersection_sum / static_cast<double>(embedding_x.size()));
//     }, threads_sizet);
//   }
//
//   // Sort the vector based on the first element (int) in ascending order
//   std::sort(IC_curve.begin(), IC_curve.end(),
//             [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
//               return a.first < b.first;
//             });
//
//   // Remove all pairs where the double value is >= 1
//   IC_curve.erase(std::remove_if(IC_curve.begin(), IC_curve.end(),
//                                 [](const std::pair<int, double>& p) {
//                                   return p.second >= 1;  // Remove if double is >= 1
//                                 }), IC_curve.end());
//
//   // Take out IC value separately
//   std::vector<double> IC_Values;
//   for (const auto& result : IC_curve) {
//     IC_Values.push_back(result.second);
//   }
//
//   // Counter to track occurrences of 1.0
//   int count = 0;
//
//   // Remove only the second occurrence of 1.0
//   auto it = std::remove_if(IC_Values.begin(), IC_Values.end(),
//                            [&count](double value) {
//                              if (value == 1.0) {
//                                count++;
//                                return count == 2;  // Remove only the second occurrence
//                              }
//                              return false;
//                            });
//
//   // Erase the marked element
//   if (it != IC_Values.end()) {
//     IC_Values.erase(it, IC_Values.end());
//   }
//
//   // Compute AUC and normalize causal strength score
//   double auc = 0.0;
//   double r_size = static_cast<double>(IC_Values.size()) - 1;
//
//   for (size_t i = 1; i < IC_Values.size(); ++i) {
//     double width = 1.0;  // Assuming uniform step size in r
//     auc += (IC_Values[i] + IC_Values[i - 1]) * width / 2.0;
//   }
//
//   double AUC = auc / r_size;
//   double ICC_score = std::max(0.0, 2.0 * (AUC - 0.5));
//
//   return ICC_score;
// }

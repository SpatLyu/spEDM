// #include <vector>
// #include <cmath>
// #include <algorithm>
// #include <numeric>
// #include <limits>
// #include <utility>
// #include <unordered_set>
//
// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score.
//  *
//  * Parameters:
//  *   - embedding_x:   The state-space reconstructed from the potential cause variable.
//  *   - embedding_y:   The state-space reconstructed from the potential effect variable.
//  *   - num_neighbors: Number of neighbors used for cross-mapping.
//  *   - n_excluded:    Number of excluded neighbors in the distance matrix.
//  *
//  * Returns:
//  *   - A double representing the CMC causal strength score, normalized between [0,1].
//  */
// double CrossMappingCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     int num_neighbors,
//     int n_excluded = 0
// ) {
//   // Ensure valid input dimensions
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return 0.0;
//   }
//
//   // Helper function to compute the distance matrix
//   auto compute_distance_matrix = [](const std::vector<std::vector<double>>& embedding) {
//     size_t n = embedding.size();
//     std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, 0.0));
//     for (size_t i = 0; i < n; ++i) {
//       for (size_t j = i + 1; j < n; ++j) {
//         double dist = 0.0;
//         for (size_t k = 0; k < embedding[i].size(); ++k) {
//           dist += std::pow(embedding[i][k] - embedding[j][k], 2);
//         }
//         distance_matrix[i][j] = distance_matrix[j][i] = std::sqrt(dist);
//       }
//     }
//     return distance_matrix;
//   };
//
//   // Compute distance matrices for embedding_x and embedding_y
//   auto distance_matrix_x = compute_distance_matrix(embedding_x);
//   auto distance_matrix_y = compute_distance_matrix(embedding_y);
//
//   // Helper function to find k-nearest neighbors
//   auto find_k_nearest_neighbors = [](const std::vector<std::vector<double>>& distance_matrix, size_t idx, int k) {
//     std::vector<size_t> neighbors(distance_matrix.size());
//     std::iota(neighbors.begin(), neighbors.end(), 0);
//     std::sort(neighbors.begin(), neighbors.end(), [&](size_t a, size_t b) {
//       return distance_matrix[idx][a] < distance_matrix[idx][b];
//     });
//     return std::vector<size_t>(neighbors.begin() + 1, neighbors.begin() + k + 1); // Exclude self
//   };
//
//   // Compute the mapping ratios
//   std::vector<double> mapping_ratios;
//   for (size_t i = 0; i < embedding_x.size(); ++i) {
//     auto neighbors_x = find_k_nearest_neighbors(distance_matrix_x, i, num_neighbors);
//     auto neighbors_y = find_k_nearest_neighbors(distance_matrix_y, i, num_neighbors + n_excluded);
//
//     // Map neighbors from embedding_x to embedding_y
//     std::unordered_set<size_t> mapped_neighbors;
//     for (size_t nx : neighbors_x) {
//       auto mapped = find_k_nearest_neighbors(distance_matrix_y, nx, 1);
//       if (!mapped.empty()) {
//         mapped_neighbors.insert(mapped[0]);
//       }
//     }
//
//     // Compute the intersection count
//     double intersection_count = 0;
//     for (size_t ny : neighbors_y) {
//       if (mapped_neighbors.find(ny) != mapped_neighbors.end()) {
//         intersection_count++;
//       }
//     }
//
//     mapping_ratios.push_back(intersection_count / neighbors_y.size());
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
//   // Normalize the score to [0,1]
//   double cmc_score = std::max(0.0, 2.0 * (auc - 0.5));
//   return cmc_score;
// }


#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "CppStats.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

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
    if (!x_all_nan || !y_all_nan) {
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

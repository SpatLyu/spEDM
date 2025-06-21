#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for grid data using simplex projection.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbors to use for prediction.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given embedding dimension.
 */
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& source,
                                              const std::vector<std::vector<double>>& target,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<int>& E,
                                              const std::vector<int>& b,
                                              int tau,
                                              int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Remove duplicates from E and b
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate all unique (E, b) pairs
  std::vector<std::pair<int, int>> unique_Ebcom;
  unique_Ebcom.reserve(Es.size() * bs.size());
  for (int e : Es) {
    for (int bn : bs) {
      unique_Ebcom.emplace_back(e, bn);
    }
  }

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(5));

  // Parallel loop over combinations
  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t i) {
    const int cur_E = unique_Ebcom[i].first;
    const int cur_b = unique_Ebcom[i].second;

    // Generate embedding
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, cur_E, tau);

    // Evaluate performance
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, cur_b);

    // Store results
    result[i][0] = cur_E;
    result[i][1] = cur_b;
    result[i][2] = metrics[0];
    result[i][3] = metrics[1];
    result[i][4] = metrics[2];
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for grid data using the S-mapping method.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& source,
                                           const std::vector<std::vector<double>>& target,
                                           const std::vector<int>& lib_indices,
                                           const std::vector<int>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E,
                                           int tau,
                                           int b,
                                           int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Generate embedding once
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, E, tau);

  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    std::vector<double> metrics = SMapBehavior(embeddings, vec_std, lib_indices, pred_indices, b, theta[i]);

    result[i][0] = theta[i];
    result[i][1] = metrics[0];
    result[i][2] = metrics[1];
    result[i][3] = metrics[2];
  }, threads_sizet);

  return result;
}

/**
 * @brief Evaluate intersection cardinality (IC) for spatial grid embeddings.
 *
 * This function computes the intersection cardinality between the k-nearest neighbors
 * of grid-embedded source and target spatial variables, across a range of embedding dimensions (E)
 * and neighborhood sizes (b). The result is an AUC (Area Under the Curve) score for each (E, b) pair
 * that quantifies the directional similarity or interaction between the spatial fields.
 *
 * The method constructs delay-like embeddings over grid cells using spatial neighborhoods,
 * filters out invalid prediction locations (e.g., with all NaN values), computes nearest neighbors
 * in embedding space, and calculates the cardinality of overlapping neighbors. These overlaps are
 * then evaluated using a CMC-based statistical test (via AUC).
 *
 * Supports both single-threaded and parallel execution using `RcppThread`.
 *
 * @param source 2D spatial variable (grid) used as the source for embedding.
 * @param target 2D spatial variable (grid) used as the target for embedding.
 * @param lib_indices Indices of spatial locations used as the library set (training).
 * @param pred_indices Indices of spatial locations used as the prediction set (evaluation).
 * @param E Vector of spatial embedding dimensions to evaluate (e.g., neighborhood sizes).
 * @param b Vector of neighbor counts (k) used to compute IC.
 * @param tau Spatial embedding spacing (lag). Determines distance between embedding neighbors.
 * @param exclude Number of nearest neighbors to exclude in IC computation.
 * @param threads Maximum number of threads to use.
 * @param parallel_level If > 0, enables parallel evaluation of b for each E.
 *
 * @return A matrix of size (|E| × |b|) × 3 with rows: [E, b, AUC]
 */
std::vector<std::vector<double>> IC4Grid(const std::vector<std::vector<double>>& source,
                                         const std::vector<std::vector<double>>& target,
                                         const std::vector<int>& lib_indices,
                                         const std::vector<int>& pred_indices,
                                         const std::vector<int>& E,
                                         const std::vector<int>& b,
                                         int tau,
                                         int exclude,
                                         int threads,
                                         int parallel_level) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Remove duplicates from E and b
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate all unique (E, b) pairs
  std::vector<std::pair<int, int>> unique_Ebcom;
  unique_Ebcom.reserve(Es.size() * bs.size());
  for (int e : Es) {
    for (int bn : bs) {
      unique_Ebcom.emplace_back(e, bn);
    }
  }

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(3));

  if (parallel_level == 0){
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = GenGridEmbeddings(source, Es[i], tau);
      auto embedding_y = GenGridEmbeddings(target, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<int> valid_pred;
      for (int idx : pred_indices) {
        if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Precompute neighbors
      auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, false, true));
      auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, false, true));

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      for (size_t j = 0; j < bs.size(); ++j){
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
          nx,ny,static_cast<size_t>(lib_indices.size()),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,0
        );

        double auc = 0;
        if (!res.empty()) auc = CppCMCTest(res[0].Intersection,">")[0];

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = auc;    // AUC
      }
    }
  } else {
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = GenGridEmbeddings(source, Es[i], tau);
      auto embedding_y = GenGridEmbeddings(target, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<int> valid_pred;
      for (int idx : pred_indices) {
        if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Precompute neighbors
      auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, false, true));
      auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, false, true));

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      RcppThread::parallelFor(0, bs.size(), [&](size_t j) {
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
          nx,ny,static_cast<size_t>(lib_indices.size()),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,1
        );

        double auc = 0;
        if (!res.empty()) auc = CppCMCTest(res[0].Intersection,">")[0];

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = auc;    // AUC
      }, threads_sizet);
    }
  }

  return result;
}

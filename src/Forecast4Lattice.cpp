#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for lattice data using simplex projection.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - nb_vec: A 2D vector of neighbor indices.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbor values to evaluate.
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given combination of E and b.
 */
std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& source,
                                                 const std::vector<double>& target,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<int>& lib_indices,
                                                 const std::vector<int>& pred_indices,
                                                 const std::vector<int>& E,
                                                 const std::vector<int>& b,
                                                 int tau,
                                                 int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Unique sorted embedding dimensions and neighbor values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate unique (E, b) combinations
  std::vector<std::pair<int, int>> unique_Ebcom;
  for (int e : Es)
    for (int bb : bs)
      unique_Ebcom.emplace_back(e, bb);

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(5));

  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t i) {
    const int Ei = unique_Ebcom[i].first;
    const int bi = unique_Ebcom[i].second;

    auto embeddings = GenLatticeEmbeddings(source, nb_vec, Ei, tau);
    auto metrics = SimplexBehavior(embeddings, target, lib_indices, pred_indices, bi);

    result[i][0] = Ei;
    result[i][1] = bi;
    result[i][2] = metrics[0]; // rho
    result[i][3] = metrics[1]; // MAE
    result[i][4] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for lattice data using the s-mapping method.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - nb_vec: A 2D vector of neighbor indices.
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
std::vector<std::vector<double>> SMap4Lattice(const std::vector<double>& source,
                                              const std::vector<double>& target,
                                              const std::vector<std::vector<int>>& nb_vec,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<double>& theta,
                                              int E,
                                              int tau,
                                              int b,
                                              int threads){
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  auto embeddings = GenLatticeEmbeddings(source, nb_vec, E, tau);
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    auto metrics = SMapBehavior(embeddings, target, lib_indices, pred_indices, b, theta[i]);

    result[i][0] = theta[i];   // theta
    result[i][1] = metrics[0]; // rho
    result[i][2] = metrics[1]; // MAE
    result[i][3] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/**
 * Compute Intersection Cardinality AUC over Lattice Embedding Settings.
 *
 * This function computes the causal strength between two lattice-structured time series
 * (`source` and `target`) by evaluating the Intersection Cardinality (IC) curve, and
 * summarizing it using the Area Under the Curve (AUC) metric.
 *
 * For each combination of embedding dimension `E` and neighbor size `b`, the function:
 *  - Generates state-space embeddings based on lattice neighborhood topology.
 *  - Filters out prediction points with missing (NaN) values.
 *  - Computes neighbor structures and evaluates intersection sizes between the mapped
 *    neighbors of `source` and `target`.
 *  - Aggregates the IC curve and estimates the AUC (optionally using significance test).
 *
 * @param source         Time series values of the potential cause variable (flattened lattice vector).
 * @param target         Time series values of the potential effect variable (same shape as `source`).
 * @param nb_vec         Neighborhood topology vector for the lattice structure.
 * @param lib_indices    Indices used for library (training) data.
 * @param pred_indices   Indices used for prediction (testing) data.
 * @param E              Vector of embedding dimensions to try.
 * @param b              Vector of neighbor sizes to try.
 * @param tau            Embedding delay (usually 1 for lattice).
 * @param exclude        Number of nearest neighbors to exclude (e.g., temporal or spatial proximity).
 * @param threads        Number of threads for parallel computation.
 * @param parallel_level Flag indicating whether to use multi-threading (0: serial, 1: parallel).
 *
 * @return A vector of size `E.size() * b.size()`, each element is a vector:
 *         [embedding_dimension, neighbor_size, auc_value].
 *         If inputs are invalid or no prediction point is valid, the AUC value is NaN.
 *
 * @note
 *   - Only AUC and p value are returned in current version. Use other utilities to derive CI.
 *   - Library and prediction indices should be adjusted for 0-based indexing before calling.
 *   - Lattice embedding assumes neighborhood-based spatial structure.
 */
std::vector<std::vector<double>> IC4Lattice(const std::vector<double>& source,
                                            const std::vector<double>& target,
                                            const std::vector<std::vector<int>>& nb_vec,
                                            const std::vector<size_t>& lib_indices,
                                            const std::vector<size_t>& pred_indices,
                                            const std::vector<int>& E,
                                            const std::vector<int>& b,
                                            int tau,
                                            int exclude,
                                            int threads,
                                            int parallel_level) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Unique sorted embedding dimensions and neighbor values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate unique (E, b) combinations
  std::vector<std::pair<int, int>> unique_Ebcom;
  for (int e : Es)
    for (int bb : bs)
      unique_Ebcom.emplace_back(e, bb);

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(4));

  if (parallel_level == 0){
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = GenLatticeEmbeddings(source, nb_vec, Es[i], tau);
      auto embedding_y = GenLatticeEmbeddings(target, nb_vec, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Precompute neighbors
      auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, false, true),lib_indices);
      auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, false, true),lib_indices);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      for (size_t j = 0; j < bs.size(); ++j){
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,0
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = cs[0];  // AUC
        result[j + bs.size() * i][3] = cs[1];  // P value
      }
    }
  } else {
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = GenLatticeEmbeddings(source, nb_vec, Es[i], tau);
      auto embedding_y = GenLatticeEmbeddings(target, nb_vec, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Precompute neighbors
      auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, false, true),lib_indices);
      auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, false, true),lib_indices);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      RcppThread::parallelFor(0, bs.size(), [&](size_t j) {
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,1
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = cs[0];  // AUC
        result[j + bs.size() * i][3] = cs[1];  // P value
      }, threads_sizet);
    }
  }

  return result;
}

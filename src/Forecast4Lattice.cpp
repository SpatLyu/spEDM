#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <tuple>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include "PatternCausality.h"
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
 *   - tau: The spatial lag step for constructing lagged state-space vectors. Default is 1.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it).  Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
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
                                                 const std::vector<int>& tau,
                                                 int style = 1,
                                                 int dist_metric = 2,
                                                 bool dist_average = true,
                                                 int threads = 8) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es)
    for (int bb : bs)
      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  RcppThread::parallelFor(0, unique_EbTau.size(), [&](size_t i) {
    const int Ei   = std::get<0>(unique_EbTau[i]);
    const int bi   = std::get<1>(unique_EbTau[i]);
    const int taui = std::get<2>(unique_EbTau[i]);
    // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

    auto embeddings = GenLatticeEmbeddings(source, nb_vec, Ei, taui, style);
    auto metrics = SimplexBehavior(embeddings, target, lib_indices, pred_indices, bi, dist_metric, dist_average);

    result[i][0] = Ei; // E
    result[i][1] = bi; // k
    result[i][2] = taui; // tau
    result[i][3] = metrics[0]; // rho
    result[i][4] = metrics[1]; // MAE
    result[i][5] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for lattice data using simplex projection (composite embeddings version).
 */
std::vector<std::vector<double>> Simplex4LatticeCom(const std::vector<double>& source,
                                                    const std::vector<double>& target,
                                                    const std::vector<std::vector<int>>& nb_vec,
                                                    const std::vector<int>& lib_indices,
                                                    const std::vector<int>& pred_indices,
                                                    const std::vector<int>& E,
                                                    const std::vector<int>& b,
                                                    const std::vector<int>& tau,
                                                    int style = 1,
                                                    int dist_metric = 2,
                                                    bool dist_average = true,
                                                    int threads = 8) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es)
    for (int bb : bs)
      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t i) {
    const int Ei   = std::get<0>(unique_EbTau[i]);
    const int bi   = std::get<1>(unique_EbTau[i]);
    const int taui = std::get<2>(unique_EbTau[i]);
    // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

    auto embeddings = GenLatticeEmbeddingsCom(source, nb_vec, Ei, taui, style);
    auto metrics = SimplexBehavior(embeddings, target, lib_indices, pred_indices, bi, dist_metric, dist_average);

    result[i][0] = Ei; // E
    result[i][1] = bi; // k
    result[i][2] = taui; // tau
    result[i][3] = metrics[0]; // rho
    result[i][4] = metrics[1]; // MAE
    result[i][5] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for lattice data using
 * the s-mapping method.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - nb_vec: A 2D vector of neighbor indices.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate. Default is 3.
 *   - tau: The spatial lag step for constructing lagged state-space vectors. Default is 1.
 *   - b: Number of nearest neighbors to use for prediction. Default is 4.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it).  Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
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
                                              int E = 3,
                                              int tau = 1,
                                              int b = 4,
                                              int style = 1,
                                              int dist_metric = 2,
                                              bool dist_average = true,
                                              int threads = 8){
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  auto embeddings = GenLatticeEmbeddings(source, nb_vec, E, tau, style);
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    auto metrics = SMapBehavior(embeddings, target, lib_indices, pred_indices, b, theta[i], dist_metric, dist_average);

    result[i][0] = theta[i];   // theta
    result[i][1] = metrics[0]; // rho
    result[i][2] = metrics[1]; // MAE
    result[i][3] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for lattice data using
 * the s-mapping method (composite embeddings version).
 */
std::vector<std::vector<double>> SMap4LatticeCom(const std::vector<double>& source,
                                                 const std::vector<double>& target,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<int>& lib_indices,
                                                 const std::vector<int>& pred_indices,
                                                 const std::vector<double>& theta,
                                                 int E = 3,
                                                 int tau = 1,
                                                 int b = 4,
                                                 int style = 1,
                                                 int dist_metric = 2,
                                                 bool dist_average = true,
                                                 int threads = 8){
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  auto embeddings = GenLatticeEmbeddingsCom(source, nb_vec, E, tau, style);
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    auto metrics = SMapBehavior(embeddings, target, lib_indices, pred_indices, b, theta[i], dist_metric, dist_average);

    result[i][0] = theta[i];   // theta
    result[i][1] = metrics[0]; // rho
    result[i][2] = metrics[1]; // MAE
    result[i][3] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

/**
 * Compute Intersection Cardinality AUC over spatial lattice data.
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
 * @param style          Embedding style selector (0: includes current state, 1: excludes it).
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads        Number of threads for parallel computation.
 * @param parallel_level Flag indicating whether to use multi-threading (0: serial, 1: parallel).
 *
 * @return A vector of size `E.size() * b.size()`, each element is a vector:
 *         [embedding_dimension, neighbor_size, auc_value, p value].
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
                                            int tau = 1,
                                            int exclude = 0,
                                            int style = 1,
                                            int dist_metric = 2,
                                            int threads = 8,
                                            int parallel_level = 0) {
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

  size_t max_num_neighbors = 0;
  if (!bs.empty()) {
    max_num_neighbors = static_cast<size_t>(bs.back() + exclude);
  }

  if (parallel_level == 0){
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = GenLatticeEmbeddings(source, nb_vec, Es[i], tau, style);
      auto embedding_y = GenLatticeEmbeddings(target, nb_vec, Es[i], tau, style);

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

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

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
      auto embedding_x = GenLatticeEmbeddings(source, nb_vec, Es[i], tau, style);
      auto embedding_y = GenLatticeEmbeddings(target, nb_vec, Es[i], tau, style);

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

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

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

/**
 * @brief Compute pattern-based causality across a parameter lattice of embedding dimension (E),
 *        neighbor count (b), and time lag (tau).
 *
 * This function builds lattice embeddings from two time series (or lattice signals) `source`
 * and `target` for every unique combination of E, b, and tau, runs pattern-based causality
 * analysis for each combination, and returns a numeric summary for all combinations.
 *
 * Steps:
 * 1. Deduplicate and sort the input parameter lists E, b, tau to obtain unique values.
 * 2. Form the Cartesian product (E, b, tau) as the search lattice.
 * 3. For each lattice node:
 *    - Generate lattice embeddings (shadow manifolds) for source and target using
 *      `GenLatticeEmbeddings(source, nb_vec, E, tau, style)`.
 *    - Call `PatternCausality(Mx, My, lib_indices, pred_indices, b, ...)` to compute
 *      pattern-level causality metrics for that node.
 * 4. Collect and return a 2D numeric matrix where each row corresponds to one lattice node
 *    and columns summarize the lattice parameters and causality metrics (e.g. TotalPos,
 *    TotalNeg, TotalDark).
 *
 * Concurrency:
 * - The function supports optional parallel execution over the lattice nodes via
 *   `RcppThread::parallelFor`. When running nodes in parallel, the inner `PatternCausality`
 *   call is forced to single-threaded to avoid nested parallelism (unless you specifically
 *   want otherwise).
 *
 * Parameters (short):
 * - source, target: raw vectors to embed.
 * - nb_vec: neighborhood structure used by the lattice embedding generator.
 * - lib_indices, pred_indices: indices used by PatternCausality for library/prediction.
 * - E, b, tau: vectors of candidate embedding dims, neighbor counts, and lags.
 * - style: embedding style passed to GenLatticeEmbeddings.
 * - zero_tolerance, dist_metric, relative, weighted: forwarded to PatternCausality.
 * - threads: hint for internal multi-threading (used in PatternCausality or parallelFor).
 * - parallel_level: 0 = low-level parallel, non-zero = high-level parallel.
 *
 * Return:
 * - matrix (std::vector<std::vector<double>>) with one row per (E,b,tau) and columns:
 *   [E, b, tau, TotalPos, TotalNeg, TotalDark]
 */
std::vector<std::vector<double>> PC4Lattice(const std::vector<double>& source,
                                            const std::vector<double>& target,
                                            const std::vector<std::vector<int>>& nb_vec,
                                            const std::vector<size_t>& lib_indices,
                                            const std::vector<size_t>& pred_indices,
                                            const std::vector<int>& E,
                                            const std::vector<int>& b,
                                            const std::vector<int>& tau,
                                            int style = 1,
                                            int zero_tolerance = 0,
                                            int dist_metric = 2,
                                            bool relative = true,
                                            bool weighted = true,
                                            int threads = 8,
                                            int parallel_level = 0) {
  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es)
    for (int bb : bs)
      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  if (parallel_level == 0){
    for (size_t i = 0; i < unique_EbTau.size(); ++i) {
      const int Ei   = std::get<0>(unique_EbTau[i]);
      const int bi   = std::get<1>(unique_EbTau[i]);
      const int taui = std::get<2>(unique_EbTau[i]);
      // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

      auto Mx = GenLatticeEmbeddings(source, nb_vec, Ei, taui, style);
      auto My = GenLatticeEmbeddings(target, nb_vec, Ei, taui, style);

      PatternCausalityRes res = PatternCausality(
        Mx, My, lib_indices, pred_indices, bi, zero_tolerance,
        dist_metric, relative, weighted, threads);

      result[i][0] = Ei;
      result[i][1] = bi;
      result[i][2] = taui;
      result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
      result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
      result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
    }
  } else {
    // Configure threads
    size_t threads_sizet = static_cast<size_t>(std::abs(threads));
    threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

    RcppThread::parallelFor(0, unique_EbTau.size(), [&](size_t i) {
      const int Ei   = std::get<0>(unique_EbTau[i]);
      const int bi   = std::get<1>(unique_EbTau[i]);
      const int taui = std::get<2>(unique_EbTau[i]);
      // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

      auto Mx = GenLatticeEmbeddings(source, nb_vec, Ei, taui, style);
      auto My = GenLatticeEmbeddings(target, nb_vec, Ei, taui, style);

      PatternCausalityRes res = PatternCausality(
        Mx, My, lib_indices, pred_indices, bi, zero_tolerance,
        dist_metric, relative, weighted, 1);

      result[i][0] = Ei;
      result[i][1] = bi;
      result[i][2] = taui;
      result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
      result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
      result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
    }, threads_sizet);
  }

  return result;
}

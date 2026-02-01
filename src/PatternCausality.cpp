#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <utility> // for std::move
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cstdint>
#include <iterator>
#include <random> // for std::mt19937_64, std::seed_seq
#include <memory> // for std::unique_ptr, std::make_unique
#include "NumericUtils.h"
#include "SymbolicDynamics.h"
#include "DataStruct.h"
#include "CppDistances.h"
#include "SignatureProjection.h"
#include <RcppThread.h>

/**
 * @brief Perform symbolic pattern–based causality analysis between real and predicted signatures.
 *
 * This function implements a deterministic, pattern–indexed causal analysis pipeline.
 * Numerical signature vectors (SMx, SMy, pred_SMy) are first transformed into symbolic
 * pattern strings using GenPatternSpace(). The function then:
 *
 *  1. Collects all unique patterns appearing in X, Y_real, and Y_pred.
 *  2. Removes patterns containing '0' (invalid placeholder state).
 *  3. Augments the pattern set by including their symmetric-opposite counterparts
 *     (swapping '1' <-> '3'), ensuring anti-diagonal causality is always representable.
 *  4. Produces a sorted, dense pattern index (0 … K-1) for deterministic and
 *     reproducible matrix alignment.
 *  5. Computes a K×K causal strength matrix M(i, j), where:
 *         i = pattern index of X(t)
 *         j = pattern index of predicted Y(t)
 *     and the matrix cell accumulates weighted causal strength over all samples.
 *
 *  6. Per-sample causality classification is performed:
 *       - No causality: pattern(Y_pred) != pattern(Y_real)
 *       - Positive     : i == j (main diagonal)
 *       - Negative     : i + j == K - 1 (anti-diagonal)
 *       - Dark         : all other off-diagonal relationships
 *
 *  7. Causal strength is optional weighted by:
 *         erf( ||pred_Y|| / (||X|| + 1e-6) )
 *     which bounds the strength in [0, 1] and prevents division instability.
 *
 *  8. A final normalized heatmap (cell-wise average) and aggregated metrics
 *     (mean positive / negative / dark strengths) are returned.
 *
 * This implemention is optimized for:
 *   - Zero-copy pattern referencing
 *   - No dynamic string comparison in the main loop
 *   - Deterministic ordering and symmetric space closure
 *   - Minimal heap reallocations
 *
 * @details
 * The key conceptual guarantee is *pattern space completeness*:
 * for any observed pattern p, the function ensures its symmetric-opposite
 * exists in the index set, even if never observed in the data. This creates
 * a fully defined anti-diagonal causal relation space, solving the correctness
 * issue in earlier implementations where i+j==K-1 could not be relied upon.
 *
 * ---------------------------------------------------------------------------
 *
 * @param SMx        X signature matrix (n × d)
 * @param SMy        Y real signature matrix (n × d)
 * @param pred_SMy   Y predicted signatures (n × d)
 * @param weighted   Whether to weight causal strength by erf(norm(pred_Y)/norm(X))
 *
 * ---------------------------------------------------------------------------
 *
 * @return PatternCausalityRes
 *
 * The result struct contains:
 *
 *   - std::vector<double> NoCausality, PositiveCausality,
 *                         NegativeCausality, DarkCausality
 *       Per-sample causal strengths (or 1 for no causality).
 *   - std::vector<int> RealLoop
 *       Indices of samples actually used (patterns valid & non-zero).
 *
 *   - std::vector<int> PatternTypes
 *       Encoded per-sample causal class:
 *         0=no causality, 1=positive, 2=negative, 3=dark.
 *
 *   - std::vector<std::string> PatternStrings
 *       Mapping index → pattern string for each row/column of the heatmap.
 *
 *   - std::vector<std::vector<double>> matrice
 *       Normalized K×K causal heatmap M(i,j).
 *
 *   - double TotalPos, TotalNeg, TotalDark
 *       Mean strength across main diagonal, anti-diagonal, and off-diagonal cells.
 *
 * ---------------------------------------------------------------------------
 *
 * @note
 *  - Patterns containing '0' are discarded from analysis.
 *  - The function guarantees symmetric pattern-space completion.
 *  - Index ordering is deterministic and reproducible across runs.
 *  - Handles NaN values robustly by ignoring them in norms and averages.
 *
 */
PatternCausalityRes GenPatternCausality(
    const std::vector<std::vector<double>>& SMx,
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& pred_SMy,
    bool weighted = true
) {
  PatternCausalityRes res;
  const size_t n = SMx.size();
  if (n == 0) return res;

  // --- 1. Generate symbolic pattern strings ---
  // GenPatternSpace() should convert numeric sequences into symbolic strings
  // (e.g., "321", "122", etc.), possibly removing or keeping NaN based on NA_rm.
  std::vector<std::string> PMx = GenPatternSpace(SMx, true);
  std::vector<std::string> PMy = GenPatternSpace(SMy, true);
  std::vector<std::string> pred_PMy = GenPatternSpace(pred_SMy, true);

  // // Basic consistency check
  // if (PMx.size() != PMy.size() || PMx.size() != pred_PMy.size()) return res;

  // --- 2. Collect unique pattern strings and mapping ---
  // Use unordered_set to collect all unique patterns efficiently
  std::unordered_set<std::string> uniq_set;
  uniq_set.reserve(PMx.size() + PMy.size() + pred_PMy.size());

  // Insert all patterns into the set
  uniq_set.insert(PMx.begin(), PMx.end());
  uniq_set.insert(PMy.begin(), PMy.end());
  uniq_set.insert(pred_PMy.begin(), pred_PMy.end());

  // Lambda to check if a pattern contains '0'
  auto str_contains_zero = [](const std::string& p) {
    return p.find('0') != std::string::npos;
  };

  // Filter out patterns containing '0' and collect them
  std::vector<std::string> filtered_patterns;
  filtered_patterns.reserve(uniq_set.size());
  for (const auto& p : uniq_set) {
    if (!str_contains_zero(p)) {
      filtered_patterns.push_back(p);
    }
  }

  // Lambda to get opposite pattern (1 <-> 3, others unchanged)
  auto get_opposite_pattern = [](const std::string& pattern) -> std::string {
    std::string opposite = pattern;
    for (char& o : opposite) {
      switch (o) {
        case '1': o = '3'; break;
        case '3': o = '1'; break;
        default: break; // Keep other characters unchanged
      }
    }
    return opposite;
  };

  // Create extended patterns by adding opposite patterns of filtered patterns with deduplication
  std::unordered_set<std::string> final_set;
  final_set.reserve(filtered_patterns.size() * 2);

  for (const auto& p : filtered_patterns) {
    final_set.insert(p); // Add original filtered patterns
    std::string opposite = get_opposite_pattern(p);
    if (opposite != p) { // Avoid adding identical opposite patterns (e.g., "222")
      final_set.insert(std::move(opposite));
    }
  }

  // Convert to vector and sort for deterministic ordering
  std::vector<std::string> unique_patterns;
  unique_patterns.reserve(final_set.size());
  unique_patterns.insert(unique_patterns.end(),
                         final_set.begin(), final_set.end());
  std::sort(unique_patterns.begin(), unique_patterns.end());

  // Build pattern_indices map
  std::unordered_map<std::string, size_t> pattern_indices;
  pattern_indices.reserve(unique_patterns.size());
  for (size_t i = 0; i < unique_patterns.size(); ++i) {
    pattern_indices.emplace(unique_patterns[i], i);
  }

  const size_t hashed_num = unique_patterns.size();
  if (hashed_num == 0) return res;
  const double midpoint = static_cast<double>(hashed_num - 1) / 2.0;

  // --- 3. Initialize result structures ---
  std::vector<std::vector<double>> heatmap_accum(
      hashed_num, std::vector<double>(hashed_num, std::numeric_limits<double>::quiet_NaN()));
  std::vector<std::vector<size_t>> count_matrix(
      hashed_num, std::vector<size_t>(hashed_num, 0));

  res.PatternStrings = std::move(unique_patterns);
  res.NoCausality.assign(n, 0.0);
  res.PositiveCausality.assign(n, 0.0);
  res.NegativeCausality.assign(n, 0.0);
  res.DarkCausality.assign(n, 0.0);
  res.PatternTypes.reserve(n);
  res.RealLoop.reserve(n);

  // --- 4. Local helper lambdas for NaN-safe math ---
  auto norm_vec_ignore_nan = [](const std::vector<double>& v) -> double {
    double sum = 0.0;
    for (double x : v) {
      if (!std::isnan(x)) sum += x * x;
    }
    return std::sqrt(sum);
  };

  auto nanmean_ignore_nan = [](const std::vector<double>& vals) -> double {
    double sum = 0.0;
    std::size_t count = 0;
    for (double v : vals) {
      if (!std::isnan(v)) {
        sum += v;
        ++count;
      }
    }
    return (count > 0)
      ? (sum / static_cast<double>(count))
      : std::numeric_limits<double>::quiet_NaN();
  };

  // --- 5. Main causality loop ---
  for (size_t t = 0; t < n; ++t) {
    // if (SMx[t].empty() || SMy[t].empty() || pred_SMy[t].empty()) continue;

    const std::string& pat_x       = PMx[t];
    const std::string& pat_y_real  = PMy[t];
    const std::string& pat_y_pred  = pred_PMy[t];

    // --- Skip invalid pattern cases ---
    if (str_contains_zero(pat_x) || str_contains_zero(pat_y_real) || str_contains_zero(pat_y_pred)) continue;

    // --- Safe index lookup ---
    auto x_it = pattern_indices.find(pat_x);
    auto y_pred_it = pattern_indices.find(pat_y_pred);
    if (x_it == pattern_indices.end() || y_pred_it == pattern_indices.end()) {
      continue;
    }

    size_t i = x_it->second;
    size_t j = y_pred_it->second;

    res.RealLoop.push_back(static_cast<int>(t));

    double strength = 0.0;
    if (pat_y_pred == pat_y_real) {
      double norm_sigx = norm_vec_ignore_nan(SMx[t]) + 1e-6;
      double ratio = norm_vec_ignore_nan(pred_SMy[t]) / norm_sigx;
      strength = weighted ? std::erf(ratio) : 1.0;
    }

    // Accumulate to heatmap
    if (std::isnan(heatmap_accum[i][j])) {
      heatmap_accum[i][j] = strength;
      count_matrix[i][j] = 1;
    } else {
      heatmap_accum[i][j] += strength;
      count_matrix[i][j] += 1;
    }

    // --- 6. Classification of per-sample causality type ---
    if (doubleNearlyEqual(strength,0.0)) {
      res.NoCausality[t] = 1.0;
      res.PatternTypes.push_back(0);
    } else {
      if (i == j && !doubleNearlyEqual(static_cast<double>(i), midpoint)) {
        res.PositiveCausality[t] = strength;
        res.PatternTypes.push_back(1);
      } else if ((i + j) == (hashed_num - 1) && !doubleNearlyEqual(static_cast<double>(i), midpoint)) {
        res.NegativeCausality[t] = strength;
        res.PatternTypes.push_back(2);
      } else {
        res.DarkCausality[t] = strength;
        res.PatternTypes.push_back(3);
      }
    }
  }

  // --- 7. Normalize heatmap by counts ---
  for (size_t i = 0; i < hashed_num; ++i) {
    for (size_t j = 0; j < hashed_num; ++j) {
      if (count_matrix[i][j] > 0) {
        heatmap_accum[i][j] = heatmap_accum[i][j] / static_cast<double>(count_matrix[i][j]);
      } else {
        heatmap_accum[i][j] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // --- 8. Compute summary metrics ---
  std::vector<double> diag_vals, anti_vals, other_vals;
  diag_vals.reserve(hashed_num);
  anti_vals.reserve(hashed_num);

  for (size_t idx = 0; idx < hashed_num; ++idx) {
    double vdiag = heatmap_accum[idx][idx];
    if (!std::isnan(vdiag)) diag_vals.push_back(vdiag);

    size_t anti_j = hashed_num - 1 - idx;
    double vanti = heatmap_accum[idx][anti_j];
    if (!std::isnan(vanti)) anti_vals.push_back(vanti);

    for (size_t j = 0; j < hashed_num; ++j) {
      if (j == idx || j == anti_j) continue;
      double v = heatmap_accum[idx][j];
      if (!std::isnan(v)) other_vals.push_back(v);
    }
  }

  res.TotalPos  = nanmean_ignore_nan(diag_vals);
  res.TotalNeg  = nanmean_ignore_nan(anti_vals);
  res.TotalDark = nanmean_ignore_nan(other_vals);
  res.matrice = std::move(heatmap_accum);

  return res;
}

/**
 * @brief Compute pattern-based causality from shadow manifolds using signature and distance-based projection.
 *
 * This function performs causality analysis between two reconstructed manifolds (`Mx`, `My`)
 * based on local neighbor projection and symbolic pattern comparison. It automates the following steps:
 *
 * 1. **Distance Computation (Dx):**
 *    - Computes pairwise distances between prediction indices (`pred_indices`)
 *      and library indices (`lib_indices`) using the chosen distance metric (`L1` or `L2`).
 *    - Parallelized with `RcppThread::parallelFor` for efficiency.
 *
 * 2. **Signature Space Generation:**
 *    - Converts the manifolds `Mx` and `My` into continuous signature spaces (`SMx`, `SMy`)
 *      via `GenSignatureSpace()`.
 *    - Supports relative embedding normalization if `relative = true`.
 *
 * 3. **Signature Projection:**
 *    - Predicts target signatures (`PredSMy`) by projecting `SMy` through local neighbors in `Dx`
 *      using `SignatureProjection()`.
 *    - Neighbors are selected by `num_neighbors`, and invalid distances (NaN) are ignored.
 *
 * 4. **Causality Computation:**
 *    - Invokes `GenPatternCausality()` to compute symbolic pattern relationships between:
 *        - real X (`SMx`), real Y (`SMy`), and predicted Y (`PredSMy`)
 *    - Produces pattern-level causality metrics, classifications, and summary matrices.
 *
 * ### Parameters
 * @param Mx             Shadow manifold for variable X (n × E)
 * @param My             Shadow manifold for variable Y (n × E)
 * @param lib_indices    Indices of library samples (used for neighbor search)
 * @param pred_indices   Indices of prediction samples (to be evaluated)
 * @param num_neighbors  Number of nearest neighbors for local projection (default = 0 → auto)
 * @param zero_tolerance Maximum number of zeros tolerated in Y signatures before truncation
 * @param dist_metric    Distance metric: 1 = L1 norm (Manhattan), 2 = L2 norm (Euclidean)
 * @param relative       Whether to normalize embedding distances relative to their local mean
 * @param weighted       Whether to weight causal strength by erf(norm(pred_Y)/norm(X))
 * @param threads        Number of threads to use (default = 1; automatically capped by hardware limit)
 *
 * ### Returns
 * @return `PatternCausalityRes` containing:
 *   - Per-pattern causality strengths (positive, negative, dark, no causality)
 *   - Causality classification summary
 *   - Heatmap-like matrix representation for downstream visualization
 *
 * ### Notes
 * - Parallelization via `RcppThread` ensures thread-safe computation of pairwise distances.
 * - Distance matrix `Dx` is asymmetric (computed only for required prediction-library pairs).
 * - This function serves as the *higher-level orchestration* combining distance, projection,
 *   and pattern causality in one pipeline.
 */
PatternCausalityRes PatternCausality(
    const std::vector<std::vector<double>>& Mx,
    const std::vector<std::vector<double>>& My,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = 0,
    int zero_tolerance = 0,
    int dist_metric = 2,
    bool relative = true,
    bool weighted = true,
    int threads = 1
){
  // Configure threads (cap at hardware concurrency)
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const size_t n_obs = Mx.size();

  // Initialize distance matrix (n_obs × n_obs) filled with NaN
  std::vector<std::vector<double>> Dx(
      n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

  // Determine distance metric: true → L1 norm, false → L2 norm
  bool L1norm = (dist_metric == 1);

  // --------------------------------------------------------------------------
  // Step 1: Compute pairwise distances between prediction and library indices
  // --------------------------------------------------------------------------
  auto compute_distance = [&](size_t p) {
    size_t pi = pred_indices[p];
    for (size_t li : lib_indices) {
      double dist = CppDistance(Mx[pi], Mx[li], L1norm, true);
      if (!std::isnan(dist)) {
        Dx[pi][li] = dist;  // assign distance; no mirroring required
      }
    }
  };

  // Parallel or serial execution depending on thread configuration
  if (threads_sizet <= 1) {
    for (size_t p = 0; p < pred_indices.size(); ++p)
      compute_distance(p);
  } else {
    RcppThread::parallelFor(0, pred_indices.size(), compute_distance, threads_sizet);
  }

  // --------------------------------------------------------------------------
  // Step 2: Generate signature spaces for Mx and My
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> SMx = GenSignatureSpace(Mx, relative);
  std::vector<std::vector<double>> SMy = GenSignatureSpace(My, relative);

  // --------------------------------------------------------------------------
  // Step 3: Predict target signatures for My using local projections
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> PredSMy = SignatureProjection(
    SMy, Dx, lib_indices, pred_indices, num_neighbors, zero_tolerance, threads_sizet);

  // --------------------------------------------------------------------------
  // Step 4: Compute pattern-based causality using symbolic pattern comparison
  // --------------------------------------------------------------------------
  PatternCausalityRes res = GenPatternCausality(SMx, SMy, PredSMy, weighted);

  return res;
}

/**
 * @brief Perform robust (bootstrapped) pattern-based causality analysis across multiple library sizes.
 *
 * This function extends `PatternCausality()` by introducing both random and systematic
 * sampling strategies for robustness evaluation. It performs repeated causality
 * estimations across different library sizes (`libsizes`) and returns results organized
 * as `[3][libsizes][boot]`:
 *
 * - Dimension 0 → metric index (0=TotalPos, 1=TotalNeg, 2=TotalDark)
 * - Dimension 1 → library size
 * - Dimension 2 → bootstrap replicate
 *
 * ### Workflow
 *
 * 1. **Distance Matrix Computation**
 *    - Computes pairwise distances between `pred_indices` and `lib_indices` once,
 *      using L1 or L2 norm (depending on `dist_metric`).
 *    - Parallelized via `RcppThread::parallelFor`.
 *    - The resulting distance matrix `Dx` is reused across all bootstraps.
 *
 * 2. **Signature Space Generation**
 *    - Builds continuous signature spaces `SMx` and `SMy` for both variables
 *      using `GenSignatureSpace()`.
 *
 * 3. **Sampling & Bootstrapping**
 *    - For each library size:
 *        - If `random_sample = true`: draw `boot` random subsets (size = L)
 *          from `lib_indices` using RNG.
 *        - If `random_sample = false`: perform deterministic slicing
 *          and **force `boot = 1`** for reproducibility.
 *
 * 4. **Causality Computation**
 *    - Projects `SMy` → `PredSMy` via `SignatureProjection()`.
 *    - Computes symbolic causality with `GenPatternCausality()`.
 *    - Extracts only the metrics `TotalPos`, `TotalNeg`, and `TotalDark`.
 *
 * 5. **Output Structure**
 *    - Returns `[3][libsizes][boot]`:
 *        - Metric index 0 → TotalPos
 *        - Metric index 1 → TotalNeg
 *        - Metric index 2 → TotalDark
 *
 * ### Parameters
 * @param Mx             Shadow manifold for variable X
 * @param My             Shadow manifold for variable Y
 * @param libsizes       Candidate library sizes
 * @param lib_indices    Indices for library samples
 * @param pred_indices   Indices for prediction samples
 * @param num_neighbors  Number of nearest neighbors for projection
 * @param boot           Number of bootstrap replicates per library size
 * @param random_sample  Whether to use random bootstrap (true) or deterministic (false)
 * @param seed           Random seed for reproducibility
 * @param zero_tolerance Max zeros allowed in signatures
 * @param dist_metric    Distance metric (1 = L1, 2 = L2)
 * @param relative       Normalize embeddings relative to local mean
 * @param weighted       Weight causality by erf(norm(pred_Y)/norm(X))
 * @param threads        Number of threads for distance/projection
 * @param parallel_level Parallelism level across boot iterations
 * @param progressbar    Whether to show progress (optional)
 *
 * @return 3D vector `[3][libsizes][boot]`
 */
std::vector<std::vector<std::vector<double>>> RobustPatternCausality(
    const std::vector<std::vector<double>>& Mx,
    const std::vector<std::vector<double>>& My,
    const std::vector<size_t>& libsizes,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = 0,
    int boot = 99,
    bool random_sample = true,
    unsigned long long seed = 42,
    int zero_tolerance = 0,
    int dist_metric = 2,
    bool relative = true,
    bool weighted = true,
    int threads = 1,
    int parallel_level = 0,
    bool progressbar = false
){
  // --------------------------------------------------------------------------
  // Step 1: Configure threads and random generators
  // --------------------------------------------------------------------------
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Enforce boot = 1 for deterministic sampling
  if (!random_sample) boot = 1;

  // Prebuild 64-bit RNG pool for reproducibility
  std::vector<std::mt19937_64> rng_pool(boot);
  for (int i = 0; i < boot; ++i) {
    std::seed_seq seq{static_cast<uint64_t>(seed), static_cast<uint64_t>(i)};
    rng_pool[i] = std::mt19937_64(seq);
  }

  // --------------------------------------------------------------------------
  // Step 2: Compute pairwise distances once
  // --------------------------------------------------------------------------
  const size_t n_obs = Mx.size();
  std::vector<std::vector<double>> Dx(
      n_obs, std::vector<double>(n_obs, std::numeric_limits<double>::quiet_NaN()));

  bool L1norm = (dist_metric == 1);

  auto compute_distance = [&](size_t p) {
    size_t pi = pred_indices[p];
    for (size_t li : lib_indices) {
      double dist = CppDistance(Mx[pi], Mx[li], L1norm, true);
      if (!std::isnan(dist)) Dx[pi][li] = dist;
    }
  };

  if (threads_sizet != 1)
    for (size_t p = 0; p < pred_indices.size(); ++p) compute_distance(p);
  else
    RcppThread::parallelFor(0, pred_indices.size(), compute_distance, threads_sizet);

  // --------------------------------------------------------------------------
  // Step 3: Generate signature spaces
  // --------------------------------------------------------------------------
  std::vector<std::vector<double>> SMx = GenSignatureSpace(Mx, relative);
  std::vector<std::vector<double>> SMy = GenSignatureSpace(My, relative);

  // --------------------------------------------------------------------------
  // Step 4: Initialize results container [3][libsizes][boot]
  // --------------------------------------------------------------------------
  const size_t n_libsizes = libsizes.size();
  std::vector<std::vector<std::vector<double>>> all_results(
      3, std::vector<std::vector<double>>(n_libsizes, std::vector<double>(boot, std::numeric_limits<double>::quiet_NaN())));

  // Optional progress bar
  std::unique_ptr<RcppThread::ProgressBar> bar;
  if (progressbar)
    bar = std::make_unique<RcppThread::ProgressBar>(n_libsizes, 1);

  // --------------------------------------------------------------------------
  // Step 5: Iterate over library sizes
  // --------------------------------------------------------------------------
  for (size_t li = 0; li < n_libsizes; ++li) {
    size_t L = libsizes[li];

    auto process_boot = [&](int b) {
      std::vector<size_t> sampled_lib, sampled_pred;

      if (random_sample) {
        std::vector<size_t> shuffled_lib = lib_indices;
        std::shuffle(shuffled_lib.begin(), shuffled_lib.end(), rng_pool[b]);
        sampled_lib.assign(shuffled_lib.begin(), shuffled_lib.begin() + L);
        // sampled_pred = sampled_lib;
      } else {
        sampled_lib.assign(lib_indices.begin(), lib_indices.begin() + L);
        // sampled_pred = sampled_lib;
      }

      std::vector<std::vector<double>> PredSMy;
      if (parallel_level == 0)
        PredSMy = SignatureProjection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, threads_sizet);
      else
        PredSMy = SignatureProjection(SMy, Dx, sampled_lib, pred_indices, num_neighbors, zero_tolerance, 1);

      PatternCausalityRes res = GenPatternCausality(SMx, SMy, PredSMy, weighted);

      all_results[0][li][b] = res.TotalPos;
      all_results[1][li][b] = res.TotalNeg;
      all_results[2][li][b] = res.TotalDark;
    };

    if (parallel_level != 0)
      RcppThread::parallelFor(0, boot, process_boot, threads_sizet);
    else
      for (int b = 0; b < boot; ++b) process_boot(b);

    if (progressbar) (*bar)++;
  }

  return all_results;
}

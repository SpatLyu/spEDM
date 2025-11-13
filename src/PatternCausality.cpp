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
#include "DataStruct.h"
#include "CppDistances.h"
#include "SignatureProjection.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the Signature Space Matrix from a State Space Matrix.
 *
 * This function transforms a state space matrix into a signature space matrix by
 * computing the differences between successive elements in each row. The transformation
 * captures dynamic patterns in state space.
 *
 * For each row in the input matrix:
 * - If `relative == true`, computes relative changes: (x[i+1] - x[i]) / x[i]
 * - If `relative == false`, computes absolute changes: x[i+1] - x[i]
 *
 * The output matrix has the same number of rows as the input, but the number of columns
 * is reduced by one (i.e., output cols = input cols - 1).
 *
 * Special handling ensures:
 * - Input validation (non-empty, at least 2 columns, numeric values)
 * - When the difference between successive states is exactly zero, the signature value is set to 0.0,
 *      indicating "no change", even in relative mode (this resolves the 0/0 undefined case for 0 → 0).
 *
 * @param mat A 2D vector representing the state space matrix.
 *            Each inner vector is a row of state coordinates.
 * @param relative If true, computes relative changes; otherwise, absolute changes.
 *                 Default is true.
 * @return A 2D vector where each row contains the signature differences of the
 *         corresponding input row. The result has dimensions [n_rows] x [n_cols - 1].
 * @throws std::invalid_argument if input is empty or has fewer than 2 columns.
 */
std::vector<std::vector<double>> GenSignatureSpace(
    const std::vector<std::vector<double>>& mat,
    bool relative = true
) {
  if (mat.empty()) {
    throw std::invalid_argument("Input matrix must not be empty.");
  }

  const size_t n_rows = mat.size();
  const size_t n_cols = mat[0].size();

  if (n_cols < 2) {
    throw std::invalid_argument("State space matrix must have at least 2 columns.");
  }

  // // Validate uniform row length
  // for (size_t i = 0; i < n_rows; ++i) {
  //   if (mat[i].size() != n_cols) {
  //     throw std::domain_error("All rows must have identical column count.");
  //   }
  // }

  const size_t out_cols = n_cols - 1;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // Pre-allocate full output matrix filled with NaN
  std::vector<std::vector<double>> result(n_rows, std::vector<double>(out_cols, nan));

  // Compute signature for each row
  for (size_t i = 0; i < n_rows; ++i) {
    const auto& row = mat[i];
    auto& out_row = result[i];

    for (size_t j = 0; j < out_cols; ++j) {
      double diff = row[j + 1] - row[j];
      // Note: NaN diff values remain NaN (meaningless pattern)
      if (!std::isnan(diff)) {
        if (doubleNearlyEqual(diff,0.0)) {
          out_row[j] = 0.0;   // no change, regardless of relative or not
        } else if (relative) {
          out_row[j] = diff / row[j];
        } else {
          out_row[j] = diff;
        }
      }
    }
  }

  return result;
}

/**
 * @brief Converts a continuous signature space matrix into a discrete string-based
 *        pattern representation for causal pattern analysis.
 *
 * This function maps each numerical signature vector (a row of the input matrix)
 * into a string of categorical symbols, encoding direction and stability of change:
 *
 *   - '0' → undefined / NaN
 *   - '1' → negative change  (value < 0)
 *   - '2' → zero change      (value == 0)
 *   - '3' → positive change  (value > 0)
 *
 *
 * @param mat   Input 2D signature matrix (n × d), real-valued, may contain NaNs.
 * @param NA_rm If true, skip rows with NaNs entirely (output "0" for those rows);
 *              if false, include NaNs as '0' symbols in their pattern strings.
 *
 * @return A vector of strings, where each string encodes one discrete pattern.
 *
 * @ Behavior controlled by `NA_rm`:
 * - **NA_rm = true (default)**:
 *   - If a row contains *any* NaN, that entire row is skipped (output is "0").
 *   - Only fully valid numeric rows are encoded into symbolic strings.
 * - **NA_rm = false**:
 *   - All rows are encoded, even those containing NaNs.
 *   - NaNs are represented as '0' in the output string.
 *
 * @ Example
 * Input:
 * ```
 * mat = [
 *   [0.1, -0.2, 0.0],
 *   [NaN, 0.3, -0.1]
 * ]
 * ```
 *
 * Output (NA_rm = true):
 * ```
 * ["312", "0"]
 * ```
 *
 * Output (NA_rm = false):
 * ```
 * ["312", "031"]
 * ```
 *
 * @Notes
 * - Each row of the returned vector corresponds to one pattern instance (time point or spatial unit).
 * - This encoding is lightweight and directly compatible with downstream
 *   pattern frequency counting, hashing, or symbolic causal inference pipelines.
 * - For large-scale symbolic modeling, string concatenation offers simplicity
 *   and transparency, trading minimal performance overhead for full interpretability.
 * - Empty input returns an empty vector.
 *
 * @ Design decisions
 * - **NaN handling**: Controlled by `NA_rm`; defaults to safe filtering mode.
 * - **Compactness**: Each row pattern stored as a single `std::string`, minimizing
 *   indexing complexity and simplifying hashing.
 * - **Performance**: Uses `reserve()` to avoid repeated allocations when building strings.
 */
std::vector<std::string> GenPatternSpace(
    const std::vector<std::vector<double>>& mat,
    bool NA_rm = true
) {
  std::vector<std::string> patterns;
  if (mat.empty()) return patterns;

  const size_t n_rows = mat.size();
  const size_t n_cols = mat[0].size();
  patterns.reserve(n_rows);

  for (size_t i = 0; i < n_rows; ++i) {
    const auto& row = mat[i];
    bool has_nan = false;
    std::string pat;
    pat.reserve(n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
      double v = row[j];
      if (std::isnan(v)) {
        has_nan = true;
        pat.push_back('0');
      } else if (doubleNearlyEqual(v,0.0)) {
        pat.push_back('2');
      } else if (v > 0.0) {
        pat.push_back('3');
      } else {
        pat.push_back('1');
      }
    }

    // When NA_rm = true and row contains NaN → replace with "0"
    if (NA_rm && has_nan) {
      patterns.emplace_back("0");
    } else {
      patterns.emplace_back(std::move(pat));
    }
  }

  return patterns;
}

/**
 * @brief Compute pattern-based causality analysis between predicted and real signatures.
 *
 * This function automatically generates symbolic patterns (PMx, PMy, pred_PMy)
 * using GenPatternSpace() from input signature matrices (SMx, SMy, pred_SMy),
 * then computes a causal strength matrix and per-sample classifications.
 *
 * Supports NaN handling through `NA_rm`:
 *  - If true, samples containing NaN are removed before pattern generation
 *    → max pattern count = 3^(E-1) + 1 (same as factorial hash style).
 *  - If false, NaN is treated as a valid symbol
 *    → max pattern count = 4^(E-1).
 *
 * Causality classification follows:
 *  - Positive causality  → pattern(Y_pred) == pattern(Y_real)
 *  - Negative causality  → symmetric opposite pattern (i + j == N - 1)
 *  - Dark causality      → other off-diagonal relationships
 *  - No causality        → zero strength
 *
 * @param SMx        X signature matrix (n × d)
 * @param SMy        Y real signature matrix (n × d)
 * @param pred_SMy   Y predicted signatures (n × d)
 * @param weighted   Whether to weight causal strength by erf(norm(pred_Y)/norm(X))
 * @param NA_rm      Whether to remove NaN samples before pattern generation
 *
 * @return PatternCausalityRes containing causal matrices, summary, and classifications.
 */
PatternCausalityRes GenPatternCausality(
    const std::vector<std::vector<double>>& SMx,
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& pred_SMy,
    bool weighted = true,
    bool NA_rm = true
) {
  PatternCausalityRes res;
  const size_t n = SMx.size();
  if (n == 0) return res;

  // --- 1. Generate symbolic pattern strings ---
  // GenPatternSpace() should convert numeric sequences into symbolic strings
  // (e.g., "0321", "1220", etc.), possibly removing or keeping NaN based on NA_rm.
  std::vector<std::string> PMx = GenPatternSpace(SMx, NA_rm);
  std::vector<std::string> PMy = GenPatternSpace(SMy, NA_rm);
  std::vector<std::string> pred_PMy = GenPatternSpace(pred_SMy, NA_rm);

  // // Basic consistency check
  // if (PMx.size() != PMy.size() || PMx.size() != pred_PMy.size()) return res;

  // --- 2. Collect unique pattern strings and mapping ---
  std::unordered_set<std::string> uniq_set;
  uniq_set.insert(PMx.begin(), PMx.end());
  uniq_set.insert(PMy.begin(), PMy.end());
  uniq_set.insert(pred_PMy.begin(), pred_PMy.end());

  std::vector<std::string> unique_patterns(uniq_set.begin(), uniq_set.end());
  res.PatternStrings = unique_patterns;

  std::unordered_map<std::string, size_t> pattern_indices;
  pattern_indices.reserve(unique_patterns.size());
  for (size_t i = 0; i < unique_patterns.size(); ++i) {
    pattern_indices.emplace(unique_patterns[i], i);
  }

  const size_t hashed_num = unique_patterns.size();
  if (hashed_num == 0) return res;

  // --- 3. Initialize result structures ---
  std::vector<std::vector<double>> heatmap_accum(
      hashed_num, std::vector<double>(hashed_num, std::numeric_limits<double>::quiet_NaN()));
  std::vector<std::vector<double>> count_matrix(
      hashed_num, std::vector<double>(hashed_num, 0.0));

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

    // --- Skip invalid pattern cases based on NA_rm ---
    if (NA_rm) {
      // In NA_rm=true mode, "0" represents an all-NaN (invalid) pattern
      if (pat_x == "0" || pat_y_real == "0" || pat_y_pred == "0") continue;
    } else {
      // In NA_rm=false mode, "000...0" (E−1 zeros) represents an all-NaN pattern
      size_t dim = SMx[0].size() > 0 ? SMx[0].size() - 1 : 1;
      std::string zero_pattern(dim, '0');
      if (pat_x == zero_pattern || pat_y_real == zero_pattern || pat_y_pred == zero_pattern) continue;
    }

    // // Skip if any pattern missing
    // if (!pattern_indices.count(pat_x) ||
    //     !pattern_indices.count(pat_y_real) ||
    //     !pattern_indices.count(pat_y_pred)) {
    //     continue;
    // }

    res.RealLoop.push_back(t);
    size_t i = pattern_indices[pat_x];
    size_t j = pattern_indices[pat_y_pred];

    double strength = 0.0;
    if (pat_y_pred == pat_y_real) {
      double norm_sigx = norm_vec_ignore_nan(SMx[t]) + 1e-6;
      double ratio = norm_vec_ignore_nan(pred_SMy[t]) / norm_sigx;
      strength = weighted ? std::erf(ratio) : 1.0;
    }

    // Accumulate to heatmap
    if (std::isnan(heatmap_accum[i][j])) {
      heatmap_accum[i][j] = strength;
      count_matrix[i][j] = 1.0;
    } else {
      heatmap_accum[i][j] += strength;
      count_matrix[i][j] += 1.0;
    }

    // --- 6. Classification of per-sample causality type ---
    if (doubleNearlyEqual(strength,0.0)) {
      res.NoCausality[t] = 1.0;
      res.PatternTypes.push_back(0);
    } else {
      double midpoint = static_cast<double>(hashed_num - 1) / 2.0;
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
  res.matrice = heatmap_accum;
  for (size_t i = 0; i < hashed_num; ++i) {
    for (size_t j = 0; j < hashed_num; ++j) {
      if (count_matrix[i][j] > 0.0) {
        res.matrice[i][j] = heatmap_accum[i][j] / count_matrix[i][j];
      } else {
        res.matrice[i][j] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // --- 8. Compute summary metrics ---
  std::vector<double> diag_vals, anti_vals, other_vals;
  diag_vals.reserve(hashed_num);
  anti_vals.reserve(hashed_num);

  for (size_t idx = 0; idx < hashed_num; ++idx) {
    double vdiag = res.matrice[idx][idx];
    if (!std::isnan(vdiag)) diag_vals.push_back(vdiag);

    size_t anti_j = hashed_num - 1 - idx;
    double vanti = res.matrice[idx][anti_j];
    if (!std::isnan(vanti)) anti_vals.push_back(vanti);

    for (size_t j = 0; j < hashed_num; ++j) {
      if (j == idx || j == anti_j) continue;
      double v = res.matrice[idx][j];
      if (!std::isnan(v)) other_vals.push_back(v);
    }
  }

  res.TotalPos  = nanmean_ignore_nan(diag_vals);
  res.TotalNeg  = nanmean_ignore_nan(anti_vals);
  res.TotalDark = nanmean_ignore_nan(other_vals);

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
 * @param NA_rm          Whether to remove NaN samples before symbolic pattern generation
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
    bool NA_rm = true,
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
  PatternCausalityRes res = GenPatternCausality(SMx, SMy, PredSMy, weighted, NA_rm);

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
 * @param NA_rm          Remove NaN samples before symbolic generation
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
    unsigned int seed = 42,
    int zero_tolerance = 0,
    int dist_metric = 2,
    bool relative = true,
    bool weighted = true,
    bool NA_rm = true,
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

      PatternCausalityRes res = GenPatternCausality(SMx, SMy, PredSMy, weighted, NA_rm);

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

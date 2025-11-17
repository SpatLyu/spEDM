#ifndef PatternCausality_H
#define PatternCausality_H

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
);

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
 */
std::vector<std::string> GenPatternSpace(
    const std::vector<std::vector<double>>& mat,
    bool NA_rm = true
);

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
 * @param sorted     Whther to sort unique pattern strings for deterministic ordering
 *
 * @return PatternCausalityRes containing causal matrices, summary, and classifications.
 */
PatternCausalityRes GenPatternCausality(
    const std::vector<std::vector<double>>& SMx,
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& pred_SMy,
    bool weighted = true,
    bool NA_rm = true,
    bool sorted = false
);

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
 * @param sorted         Whther to sort unique pattern strings for deterministic ordering
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
    bool sorted = false,
    int threads = 1
);

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
);

#endif // PatternCausality_H

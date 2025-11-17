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

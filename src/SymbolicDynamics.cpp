#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include "NumericUtils.h"

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

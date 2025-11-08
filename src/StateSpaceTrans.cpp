#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <cstdint> // for uint8_t

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
        if (diff == 0.0) {
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
 * @brief Transforms a signature space matrix into a discrete pattern space matrix
 *        for causal pattern analysis and symbolic dynamics encoding.
 *
 * This function converts each real-valued element in the input signature matrix
 * into a categorical symbol based on its sign and magnitude:
 *   - 0: undefined pattern (input was NaN or mathematically indeterminate)
 *   - 1: negative change (value < 0)  → "decrease"
 *   - 2: zero change (value == 0)     → "no-change"
 *   - 3: positive change (value > 0)  → "increase"
 *
 * The mapping is designed to support downstream causal inference workflows
 * (e.g., pattern causality heatmaps, transition counting, and symbolic dynamics),
 * where continuous signature values are abstracted into interpretable discrete states.
 *
 * Key design choices:
 *   - **NaN handling**: Input NaN values (e.g., from 0/0 in relative mode) are mapped to 0,
 *     providing a consistent "invalid/undefined" marker that can be filtered out later.
 *   - **Exact zero detection**: Only values exactly equal to 0.0 are mapped to "no-change" (2).
 *     This assumes that the input signature matrix has been preprocessed such that
 *     true "no-change" states are represented as exact zeros (e.g., via diff == 0 logic
 *     in GenSignatureSpace). Floating-point noise should be handled upstream if needed.
 *   - **Memory efficiency**: Uses std::uint8_t (1 byte per element) to minimize memory footprint,
 *     as only 4 distinct states (0–3) need to be represented. This reduces memory usage
 *     by 75% compared to int32_t and 87.5% compared to size_t on 64-bit systems.
 *   - **Type safety**: Avoids floating-point types for categorical data, preventing
 *     accidental arithmetic on pattern codes and improving code clarity.
 *
 * The output matrix has the exact same dimensions as the input:
 *   - Number of rows: preserved (each row corresponds to a time point or trajectory)
 *   - Number of columns: preserved (each column corresponds to a lagged difference or embedding dimension)
 *
 * This function is typically used after GenSignatureSpace and before:
 *   - Pattern hashing (e.g., via string conversion or rolling hash)
 *   - Transition matrix construction
 *   - Causal pattern matching (e.g., in analyze_pc_causality-style algorithms)
 *
 * @param mat A 2D matrix of signature values (output from GenSignatureSpace).
 *            Expected to be a dense matrix of doubles, possibly containing NaNs.
 *            Must be non-empty and rectangular (all rows same length).
 *
 * @return A 2D matrix of type std::uint8_t with the same shape as input,
 *         where each element is an integer in {0, 1, 2, 3} representing the
 *         discrete pattern state as defined above.
 *
 * @note Empty input returns empty output (no exception).
 *
 * @see GenSignatureSpace
 */
std::vector<std::vector<std::uint8_t>> GenPatternSpace(
    const std::vector<std::vector<double>>& mat
) {
  if (mat.empty()) return {};

  const size_t n_rows = mat.size();
  const size_t n_cols = mat[0].size();

  // Preallocate full matrix with zeros
  std::vector<std::vector<std::uint8_t>> result(
      n_rows, std::vector<std::uint8_t>(n_cols, 0));

  for (size_t i = 0; i < n_rows; ++i) {
    const auto& row = mat[i];
    auto& out_row = result[i];  // direct reference to avoid multiple lookups

    for (size_t j = 0; j < n_cols; ++j) {
      double v = row[j];
      // Note: NaN values remain 0 (undefined pattern)
      if (!std::isnan(v)){
        if (v < 0.0) {
          out_row[j] = 1;   // negative change
        } else if (v > 0.0) {
          out_row[j] = 3;   // positive change
        } else {
          out_row[j] = 2;   // no change
        }
      }
    }
  }

  return result;
}

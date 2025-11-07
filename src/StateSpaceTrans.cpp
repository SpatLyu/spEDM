#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

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
 *
 * @param mat A 2D vector representing the state space matrix.
 *            Each inner vector is a row of state coordinates.
 * @param relative If true, computes relative changes; otherwise, absolute changes.
 *                 Default is true.
 * @return A 2D vector where each row contains the signature differences of the
 *         corresponding input row. The result has dimensions [n_rows] x [n_cols - 1].
 * @throws std::invalid_argument if input is empty or has fewer than 2 columns.
 */
std::vector<std::vector<double>> signatureSpace(
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
      if (relative) {
        double denom = row[j];
        if (denom != 0) { // temporally implement
          out_row[j] = diff / denom;
        }
      } else {
        out_row[j] = diff;
      }
    }
  }

  return result;
}

#ifndef SymbolicDynamics_H
#define SymbolicDynamics_H

#include <vector>
#include <cmath>
#include <limits>
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

#endif // SymbolicDynamics_H

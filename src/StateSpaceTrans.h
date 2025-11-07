#ifndef StateSpaceTrans_H
#define StateSpaceTrans_H

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
);

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
);

#endif // StateSpaceTrans_H

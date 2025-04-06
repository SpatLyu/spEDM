#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <utility>
#include "CppStats.h"

/**
 * Converts a 2D grid position (row, column) to a 1D index in row-major order.
 * This function calculates the corresponding one-dimensional index based on the given
 * row and column position within a grid of specified dimensions.
 *
 * Parameters:
 *   curRow   - The current row number (1-based indexing).
 *   curCol   - The current column number (1-based indexing).
 *   totalRow - The total number of rows in the grid.
 *   totalCol - The total number of columns in the grid.
 *
 * Returns:
 *   The computed 1D index (0-based indexing) corresponding to the input row and column.
 */
int LocateGridIndices(int curRow, int curCol,
                      int totalRow, int totalCol);

/**
 * Converts a 1D grid index (cell number) to a 2D grid position (row, column) in row-major order.
 * This function determines the corresponding row and column indices based on a given
 * one-dimensional index within a grid of specified dimensions.
 *
 * Parameters:
 *   cellNum  - The 1D index of the cell (0-based indexing).
 *   totalCol - The total number of columns in the grid.
 *
 * Returns:
 *   A std::vector<int> containing the row index and column index (both 0-based indexing).
 */
std::vector<int> RowColFromGrid(int cellNum, int totalCol);

/**
 * Converts a 2D grid data matrix (vector of vectors) to a 1D vector by concatenating the rows.
 * This function iterates over each row of the input matrix and appends the elements to the resulting vector.
 *
 * Parameters:
 *   Matrix - A 2D vector containing the grid data, where each element is a row of double values.
 *
 * Returns:
 *   A 1D vector containing all the elements of the input matrix, arranged row by row.
 */
std::vector<double> GridMat2Vec(const std::vector<std::vector<double>>& Matrix);

/**
 * Converts a 1D vector to a 2D grid data matrix by filling the matrix row by row.
 * This function assumes the total number of elements in the vector is exactly divisible by the specified number of rows.
 *
 * Parameters:
 *   Vec   - A 1D vector containing the grid data, where elements are arranged in row-major order.
 *   NROW  - The number of rows in the resulting matrix.
 *
 * Returns:
 *   A 2D vector (matrix) containing the grid data, arranged by rows.
 */
std::vector<std::vector<double>> GridVec2Mat(const std::vector<double>& Vec,
                                             int NROW);

/**
 * Computes the lagged values for each element in a grid matrix based on a specified lag number and Moore neighborhood.
 * For each element in the matrix, the function calculates the values of its neighbors at a specified lag distance
 * in each of the 8 directions of the Moore neighborhood. If a neighbor is out of bounds, it is assigned a NaN value.
 *
 * Parameters:
 *   mat    - A 2D vector representing the grid data.
 *   lagNum - The number of steps to lag when considering the neighbors in the Moore neighborhood.
 *
 * Returns:
 *   A 2D vector containing the lagged values for each element in the grid, arranged by the specified lag number.
 *   If a neighbor is out of bounds, it is filled with NaN.
 *
 * Note:
 *   The return value for each element is the lagged value of the neighbors, not the index of the neighbor.
 */
std::vector<std::vector<double>> CppLaggedVar4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
);

/**
 * Generates grid embeddings by calculating lagged variables for each element in a grid matrix,
 * and stores the results in a matrix where each row represents an element and each column represents
 * a different lagged value or the original element.
 *
 * Parameters:
 *   mat  - A 2D vector representing the grid data.
 *   E    - The number of embedding dimensions (columns in the resulting matrix).
 *   tau  - The spatial lag step for constructing lagged state-space vectors.
 *
 * Returns:
 *   A 2D vector (matrix) where each row contains the averaged lagged variables for
 *   each embedding dimension (column). Columns where all values are NaN are removed.
 *
 * Note:
 *   When tau = 0, the lagged variables are calculated for lag steps of 0, 1, ..., E-1.
 *   When tau > 0, the lagged variables are calculated for lag steps of tau, 2*tau, ..., E*tau,
 *   and this means the actual lag steps form an arithmetic sequence with a common difference of tau.
 */
std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    int tau
);

/**
 * @brief Perform grid-based symbolization on a 2D numeric matrix using local neighborhood statistics.
 *
 * This function calculates a symbolic representation (`fs`) for each non-NaN grid cell based on its
 * k most similar neighbors in terms of value difference. The process expands radially (Queen's case)
 * around each cell until at least k valid neighbors are collected. The final symbol for each cell
 * reflects its local homogeneity pattern relative to the global median.
 *
 * The steps include:
 * 1. Flatten the matrix to compute the global median `s_me` using all valid values.
 * 2. For each grid cell:
 *    - Find up to k nearest neighbors by increasing the neighborhood lag until k valid neighbors are found.
 *    - Sort neighbors by absolute difference in value from the center cell.
 *    - Select the top k values and compute a first indicator vector (`tau_s`) by comparing to global median.
 *    - Compute a second indicator vector (`l_s`) by comparing `tau_s[i]` to the center cell’s own indicator (`taus`).
 *    - Sum `l_s` to get a symbolic value `fs` representing the symbolic spatial consistency.
 *
 * @param mat A 2D grid (matrix) of values to be symbolized. `NaN` values are treated as missing.
 * @param k The number of neighbors to consider for the symbolization of each cell.
 *
 * @return A flattened 1D vector representing the symbolic values of each grid cell (row-major order).
 *         Cells with no valid value or insufficient neighbors remain as `NaN`.
 *
 * Note:
 * - Uses Queen's neighborhood definition for expanding neighborhoods (8 directions per layer).
 * - Grid edges and missing values are handled robustly during expansion.
 * - Useful for symbolic dynamics, pattern analysis, or spatial entropy estimation.
 */
std::vector<double> GenGridSymbolization(
    const std::vector<std::vector<double>>& mat,
    size_t k);

#endif // CppGridUtils_H

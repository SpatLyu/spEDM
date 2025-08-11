#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <set>
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
std::vector<std::vector<double>> CppLaggedVal4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
);

/**
 * Generates grid embeddings by calculating lagged variables for each element in a grid matrix,
 * and stores the results in a matrix where each row represents an element and each column represents
 * a different lagged value or the original element.
 *
 * Parameters:
 *   mat   - A 2D vector representing the grid data.
 *   E     - The number of embedding dimensions (columns in the resulting matrix).
 *   tau   - The spatial lag step for constructing lagged state-space vectors.
 *   style - Embedding style selector:
 *             - style = 0: embedding includes current state as the first dimension.
 *             - style != 0: embedding excludes current state.
 *
 * Returns:
 *   A 2D vector (matrix) where each row contains the averaged lagged variables for
 *   each embedding dimension (column). Columns where all values are NaN are removed.
 *
 * Note:
 *   When tau = 0, lagged variables are calculated for lag steps 0, 1, ..., E-1.
 *   When tau > 0 and style = 0, lagged variables are calculated for lag steps 0, tau, 2*tau, ..., (E-1)*tau.
 *   When tau > 0 and style != 0, lagged variables are calculated for lag steps tau, 2*tau, ..., E*tau.
 */
std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    int tau,
    int style = 1
);

/**
 * @brief Generate k nearest neighbors for all cells in a grid,
 *        choosing only from cells listed in lib and using Queen adjacency.
 *
 * This function identifies, for each cell in a 2D grid (represented by `mat`), its k nearest neighbors
 * among a set of valid cells specified by `lib`. Neighbor relationships are determined based on
 * Queen contiguity (8 directions: N, NE, E, SE, S, SW, W, NW), recursively expanding outward if fewer
 * than k valid neighbors are found in immediate adjacency. The "distance" between cells is calculated
 * as the absolute difference in cell values. When more than k valid candidates are available, the k
 * closest cells (in value) are selected.
 *
 * The output is a vector of size mat.size() × k, where each entry corresponds to the linear index
 * (row-major) of one neighbor per cell, and the output is grouped per row cell.
 *
 * @param mat 2D matrix of values (with possible NaNs).
 * @param lib List of valid cells (flattened indices) to choose neighbors from.
 * @param k   Number of neighbors to find for each cell.
 * @return    A list of neighbors for each cell in mat, in row-major order.
 */
std::vector<std::vector<int>> GenGridNeighbors(
    const std::vector<std::vector<double>>& mat,
    const std::vector<int>& lib,
    size_t k);

/**
 * @brief Perform grid-based symbolization on a 2D numeric matrix using local neighborhood statistics.
 *
 * This function computes a symbolic representation (`fs`) for each non-NaN grid cell in `pred` based on its
 * k most similar neighbors selected from the `lib` set. Each symbol reflects local spatial consistency with
 * respect to a global threshold derived from known values.
 *
 * The updated process includes:
 * 1. Compute the global threshold (`s_me`) as the median of all valid (non-NaN) values from `lib` locations only.
 * 2. For each prediction location in `pred`:
 *    - Identify neighbors from `lib` by radially expanding the neighborhood (Queen's case) until at least `k` valid neighbors are found.
 *    - Rank neighbors by absolute difference from the center cell's value.
 *    - Select the top `k` neighbors and assign each a binary symbol (`tau_s`) based on comparison to the global median (`s_me`).
 *    - Compare each `tau_s` to the center cell's own binary symbol (`taus`), and count how many are consistent.
 *    - Store this count as the final symbolic value `fs` for the cell.
 *
 * @param mat A 2D grid (matrix) of values to be symbolized. `NaN` values are treated as missing.
 * @param lib A set of valid reference (library) cell locations as (row, col) pairs. Used for both thresholding and neighbor selection.
 * @param pred A set of prediction cell locations as (row, col) pairs where symbols will be computed.
 * @param k The number of nearest neighbors (by value similarity) to use in symbolization.
 *
 * @return A 1D vector of symbolic values corresponding to the prediction locations, in row-major order.
 *         Entries remain `NaN` if the center cell is invalid or lacks enough valid neighbors.
 *
 * Notes:
 * - Uses Queen’s case (8-directional) neighborhood expansion for neighbor search.
 * - Median is computed only from valid values in the `lib` set (not the whole matrix).
 * - Particularly useful in spatial symbolic analysis, spatial entropy estimation, or causal lattice models.
 */
std::vector<double> GenGridSymbolization(
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::pair<int, int>>& lib,
    const std::vector<std::pair<int, int>>& pred,
    size_t k);

/**
 * @brief Divide a 2D grid (matrix) into approximately square or shaped blocks.
 *
 * This function partitions a 2D matrix into `b` blocks of roughly equal size.
 * The matrix is represented as a vector of row vectors and assumed to be
 * row-major (i.e., each inner vector represents a row of the matrix).
 *
 * The grid can be divided in different ways based on the `shape` parameter:
 * - 1: Horizontal cuts (blocks are divided row-wise).
 * - 2: Vertical cuts (blocks are divided column-wise).
 * - 3: Radial cuts (blocks are divided in a radial pattern from the center).
 *
 * The grid is divided by first estimating a grid layout of `br` rows and `bc` columns
 * such that br * bc >= b and the blocks are as square as possible for shapes 1 and 2.
 * For radial cuts, the grid is divided into concentric rings and sectors.
 *
 * Each cell in the matrix is assigned a block ID ranging from 0 to b-1, stored in
 * a 1D vector corresponding to the flattened row-major order of the matrix.
 *
 * Any leftover blocks are merged into the last block.
 *
 * @param mat A 2D grid represented as a vector of vectors (row-major).
 * @param b   Number of blocks to divide the grid into.
 * @param shape The shape of the cuts (1: horizontal, 2: vertical, 3: radial).
 * @return A vector of size rows * cols where each element is the block label
 *         assigned to the corresponding cell.
 */
std::vector<int> CppDivideGrid(
    const std::vector<std::vector<double>>& mat,
    int b, int shape = 3);

#endif // CppGridUtils_H

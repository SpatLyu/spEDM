#ifndef CppGridUtils_H
#define CppGridUtils_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

/**
 * Calculate the one-dimensional index for a specified position in a 2D grid.
 *
 * This function converts a row and column position in a 2D grid into a corresponding
 * one-dimensional index, assuming the grid is stored in row-major order (i.e., rows are stored sequentially).
 *
 * @param curRow   The current row number (1-based indexing).
 * @param curCol   The current column number (1-based indexing).
 * @param totalRow The total number of rows in the grid.
 * @param totalCol The total number of columns in the grid.
 * @return         The calculated one-dimensional index (0-based indexing).
 */
int LocateGridIndices(int curRow, int curCol,
                      int totalRow, int totalCol);

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
 * Generates embeddings for a given vector and neighborhood matrix by computing the mean of neighbor values
 * for each spatial unit, considering both the immediate neighbors and neighbors up to a specified lag number.
 *
 * Parameters:
 *   vec  - A vector of values, one for each spatial unit, to be embedded.
 *   nb   - A 2D matrix where each row represents the neighbors of a spatial unit.
 *   E    - The embedding dimension, specifying how many lags to consider in the embeddings.
 *   tau  - The spatial lag step for constructing lagged state-space vectors.
 *
 * Returns:
 *   A 2D vector representing the embeddings for each spatial unit, where each spatial unit has a row in the matrix
 *   corresponding to its embedding values for each lag number. If no valid embedding columns remain after removing
 *   columns containing only NaN values, a filtered matrix is returned.
 */
std::vector<std::vector<double>> GenLatticeEmbeddings(
    const std::vector<double>& vec,
    const std::vector<std::vector<int>>& nb,
    int E,
    int tau);

#endif // CppGridUtils_H

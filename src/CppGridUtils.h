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

std::vector<std::vector<double>> CppLaggedVar4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
);

std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    bool includeself);

#endif // CppGridUtils_H

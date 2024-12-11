#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

/**
 * Expands a matrix by adding rows and columns of NaN values based on the lag number.
 *
 * @param mat The input matrix represented as a 2D vector of doubles.
 * @param lagNum The number of times the matrix should be expanded. Must be a non-negative integer.
 * @return The expanded matrix.
 */
std::vector<std::vector<double>> ExpandMat(std::vector<std::vector<double>> mat,
                                           int lagNum) {
  // If lagNum is negative, return the original matrix as no expansion is needed.
  if (lagNum < 0) {
    return mat;
  }

  // If lagNum is greater than 1, recursively call ExpandMat to expand the matrix.
  if (lagNum > 1) {
    mat = ExpandMat(mat, lagNum - 1);
  }

  // Get the number of columns and rows in the matrix.
  int numCols = mat[0].size();
  // int numRows = mat.size();

  // Add a row of NaN values at the top and bottom of the matrix.
  std::vector<double> nanRow(numCols, std::numeric_limits<double>::quiet_NaN());
  mat.insert(mat.begin(), nanRow); // Add at the top
  mat.push_back(nanRow);          // Add at the bottom

  // Add a column of NaN values at the left and right of the matrix.
  for (auto& row : mat) {
    row.insert(row.begin(), std::numeric_limits<double>::quiet_NaN()); // Add at the left
    row.push_back(std::numeric_limits<double>::quiet_NaN());          // Add at the right
  }

  return mat;
}

/**
 * Generates a 2D matrix of lagged variables based on the input matrix and lag number.
 *
 * @param mat The input matrix represented as a 2D vector of doubles.
 * @param lagNum The lag number used to generate the lagged variables.
 * @return A 2D matrix of lagged variables.
 */
std::vector<std::vector<double>> LaggedVariableAs2Dim(
    std::vector<std::vector<double>> mat,
    int lagNum) {
  // Get the number of columns and rows in the input matrix
  int numCols = mat[0].size();
  int numRows = mat.size();

  // Expand the matrix using the ExpandMat function
  mat = ExpandMat(mat, lagNum);

  // Initialize the lagged variable matrix with dimensions (numRows * numCols) x (8 * lagNum)
  int laggedVarRows = numRows * numCols;
  int laggedVarCols = 8 * lagNum;
  std::vector<std::vector<double>> laggedVar(laggedVarRows, std::vector<double>(laggedVarCols, std::numeric_limits<double>::quiet_NaN()));

  // Iterate over each row and column of the original matrix
  for (int r = 1; r <= numRows; ++r) {
    for (int c = 1; c <= numCols; ++c) {
      int item = 0; // Index for the laggedVar matrix

      // Generate the sequence of column offsets
      int colsAddStart = -std::trunc((3 + (lagNum - 1) * 2) / 2);
      int colsAddEnd = std::trunc((3 + (lagNum - 1) * 2) / 2);
      for (int la = colsAddStart; la <= colsAddEnd; ++la) {
        laggedVar[(r - 1) * numCols + c - 1][item] = mat[r - 1 + lagNum][c - 1 + la + lagNum];
        item++;
      }

      // Generate the sequence of row offsets
      for (int ra = -(lagNum - 1); ra <= (lagNum - 1); ++ra) {
        laggedVar[(r - 1) * numCols + c - 1][item] = mat[r - 1 + ra + lagNum][c - 1 - lagNum + lagNum];
        item++;

        laggedVar[(r - 1) * numCols + c - 1][item] = mat[r - 1 + ra + lagNum][c - 1 + lagNum + lagNum];
        item++;
      }

      // Generate the sequence of column offsets again
      for (int la = colsAddStart; la <= colsAddEnd; ++la) {
        laggedVar[(r - 1) * numCols + c - 1][item] = mat[r - 1 + lagNum + lagNum][c - 1 + la + lagNum];
        item++;
      }
    }
  }

  return laggedVar;
}

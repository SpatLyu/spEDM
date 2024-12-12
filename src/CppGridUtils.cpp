#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

// Note that the return value is the value of the lagged order position, not the index.
std::vector<std::vector<double>> CppLaggedVar4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum
) {
  int numRows = mat.size();
  int numCols = mat[0].size();
  int totalElements = numRows * numCols;
  std::vector<std::vector<double>> result(totalElements, std::vector<double>(8, std::numeric_limits<double>::quiet_NaN()));

  // Directions for Moore neighborhood (8 directions)
  const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};
  const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      int elementIndex = r * numCols + c;
      int neighborIndex = 0;

      // Calculate positions that are exactly 'lagNum' units away
      for (int dir = 0; dir < 8; ++dir) {
        int dr = dx[dir] * lagNum;
        int dc = dy[dir] * lagNum;
        int rNeighbor = r + dr;
        int cNeighbor = c + dc;

        // Check if the neighbor position is within bounds
        if (rNeighbor >= 0 && rNeighbor < numRows && cNeighbor >= 0 && cNeighbor < numCols) {
          result[elementIndex][neighborIndex] = mat[rNeighbor][cNeighbor];
        }
        // Else, it remains NaN

        neighborIndex++;
        if (neighborIndex >= 8) {
          break;
        }
      }
    }
  }

  return result;
}

std::vector<std::vector<std::vector<double>>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat, int E) {
  // Initialize a vector to store the embeddings
  std::vector<std::vector<std::vector<double>>> xEmbeddings(E + 1);

  // The first embedding is the transpose of the input matrix
  int numRows = mat.size();
  int numCols = mat[0].size();
  xEmbeddings[0].resize(numCols, std::vector<double>(numRows));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      xEmbeddings[0][c][r] = mat[r][c]; // Transpose the matrix
    }
  }

  // Generate the remaining embeddings using laggedVariableAs2Dim
  for (int i = 1; i <= E; ++i) {
    xEmbeddings[i] = CppLaggedVar4Grid(mat, i);
  }

  return xEmbeddings;
}

std::vector<std::vector<double>> GenGridEmbeddings2(
    const std::vector<std::vector<double>>& mat,
    int E) {
  // Calculate the total number of elements in all subsets of mat
  int total_elements = 0;
  for (const auto& subset : mat) {
    total_elements += subset.size();
  }

  // Initialize the result matrix with total_elements rows and E+1 columns
  std::vector<std::vector<double>> result(total_elements, std::vector<double>(E + 1, 0.0));

  // Fill the first column with the elements from mat
  int row = 0;
  for (const auto& subset : mat) {
    for (double value : subset) {
      result[row][0] = value;
      ++row;
    }
  }

  // Fill the remaining columns (2 to E+1) with the averaged lagged variables
  for (int lagNum = 1; lagNum <= E; ++lagNum) {
    // Calculate the lagged variables for the current lagNum
    std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);

    // Fill the current column (lagNum + 1) with the averaged lagged variables
    row = 0;
    for (const auto& subset : lagged_vars) {
      double sum = 0.0;
      int count = 0;
      for (int i = 0; i < 8; ++i) {
        double val = subset[i];
        if (!std::isnan(val)) {
          sum += val;
          ++count;
        }
      }
      if (count > 0) {
        result[row][lagNum] = sum / count;
      } else {
        result[row][lagNum] = std::numeric_limits<double>::quiet_NaN();
      }
      ++row;
    }
  }

  return result;
}

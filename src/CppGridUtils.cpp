#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

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
int LocateGridIndices(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
}

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
std::vector<int> RowColFromGrid(int cellNum, int totalCol) {
  int row = cellNum / totalCol;
  int col = cellNum % totalCol;
  return {row, col};
}

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
std::vector<double> GridMat2Vec(const std::vector<std::vector<double>>& Matrix){
  std::vector<double> vec;
  for (const auto& row : Matrix) {
    vec.insert(vec.end(), row.begin(), row.end());
  }
  return vec;
}

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
                                             int NROW){
  // Calculate the number of columns based on the vector size and number of rows
  int NCOL = Vec.size() / NROW;

  // Create the resulting matrix with NROW rows and NCOL columns
  std::vector<std::vector<double>> matrix(NROW, std::vector<double>(NCOL));

  // Fill the matrix with values from the input vector
  for (int i = 0; i < NROW; ++i) {
    for (int j = 0; j < NCOL; ++j) {
      matrix[i][j] = Vec[i * NCOL + j];
    }
  }

  return matrix;
}

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
) {
  // Validate input
  if (mat.empty() || mat[0].empty() || lagNum < 0) {
    return {};
  }

  const int rows = mat.size();
  const int cols = mat[0].size();
  const int numCells = rows * cols;
  const int numNeighbors = 8 * lagNum;

  // If lagNum is 0, return the current values of gird row by row
  if (lagNum == 0) {
    std::vector<std::vector<double>> result;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        result.push_back({mat[i][j]});
      }
    }
    return result;
  }

  // Generate all valid offsets for the given lagNum (Queen's case)
  std::vector<std::pair<int, int>> offsets;
  for (int dx = -lagNum; dx <= lagNum; ++dx) {
    for (int dy = -lagNum; dy <= lagNum; ++dy) {
      if (std::max(std::abs(dx), std::abs(dy)) == lagNum) {
        offsets.emplace_back(dx, dy);
      }
    }
  }

  // Initialize result with NaN
  std::vector<std::vector<double>> result(
      numCells,
      std::vector<double>(numNeighbors, std::numeric_limits<double>::quiet_NaN())
  );

  // Populate neighbor values
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      const int cellIndex = i * cols + j;
      for (size_t k = 0; k < offsets.size(); ++k) {
        const auto& [dx, dy] = offsets[k];
        const int ni = i + dx;
        const int nj = j + dy;
        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
          result[cellIndex][k] = mat[ni][nj];
        }
        // Else remains NaN
      }
    }
  }

  return result;
}

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
) {
  int numRows = mat.size();
  int numCols = mat[0].size();
  int total_elements = numRows * numCols;

  // Initialize the result matrix with total_elements rows and E columns and fill with NaN values
  std::vector<std::vector<double>> result(total_elements, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));

  if (tau == 0) {
    // Flatten the matrix (mat) into the first column of the result matrix
    int row = 0;
    for (const auto& subset : mat) {
      for (double value : subset) {
        result[row][0] = value;
        ++row;
      }
    }

    // Fill the remaining columns (2 to E) with the averaged lagged variables
    for (int lagNum = 1; lagNum < E; ++lagNum) {
      // Calculate the lagged variables for the current lagNum
      std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);

      // Check if all elements in lagged_vars are NaN
      bool allNaN = true;
      for (const auto& subset : lagged_vars) {
        for (double val : subset) {
          if (!std::isnan(val)) {
            allNaN = false;
            break;
          }
        }
        if (!allNaN) break;
      }

      // If all elements are NaN, stop further processing for this lagNum
      if (allNaN) {
        break;
      }

      // Fill the current column (lagNum) with the averaged lagged variables
      row = 0;
      for (const auto& subset : lagged_vars) {
        double sum = 0.0;
        int count = 0;
        for (double val : subset) {
          if (!std::isnan(val)) {
            sum += val;
            ++count;
          }
        }

        if (count > 0) {
          result[row][lagNum] = sum / count; // Average the valid values
        }
        ++row;
      }
    }
  } else {
    // When tau != 0, calculate lagged variables for tau, 2*tau, ..., E*tau
    int row = 0;
    for (int i = 1; i <= E; ++i) {
      int lagNum = i * tau;  // Calculate the actual lag step

      // Calculate the lagged variables for the current lagNum
      std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);

      // Check if all elements in lagged_vars are NaN
      bool allNaN = true;
      for (const auto& subset : lagged_vars) {
        for (double val : subset) {
          if (!std::isnan(val)) {
            allNaN = false;
            break;
          }
        }
        if (!allNaN) break;
      }

      // If all elements are NaN, stop further processing for this lagNum
      if (allNaN) {
        break;
      }

      // Fill the current column (i-1) with the averaged lagged variables
      row = 0;
      for (const auto& subset : lagged_vars) {
        double sum = 0.0;
        int count = 0;
        for (double val : subset) {
          if (!std::isnan(val)) {
            sum += val;
            ++count;
          }
        }

        if (count > 0) {
          result[row][i - 1] = sum / count; // Average the valid values
        }
        ++row;
      }
    }
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns; // To store indices of valid columns

  // Iterate over each column to check if it contains any non-NaN values
  for (size_t col = 0; col < result[0].size(); ++col) {
    bool isAllNaN = true;
    for (size_t row = 0; row < result.size(); ++row) {
      if (!std::isnan(result[row][col])) {
        isAllNaN = false;
        break;
      }
    }
    if (!isAllNaN) {
      validColumns.push_back(col); // Store the index of valid columns
    }
  }

  // If no columns are removed, return the original result
  if (validColumns.size() == result[0].size()) {
    return result;
  } else {
    // // Issue a warning if any columns are removed
    // std::cerr << "Warning: remove all-NA embedding vector columns caused by excessive embedding dimension E selection." << std::endl;

    // Construct the filtered embeddings matrix
    std::vector<std::vector<double>> filteredEmbeddings;
    for (size_t row = 0; row < result.size(); ++row) {
      std::vector<double> filteredRow;
      for (size_t col : validColumns) {
        filteredRow.push_back(result[row][col]);
      }
      filteredEmbeddings.push_back(filteredRow);
    }

    // Return the filtered embeddings matrix
    return filteredEmbeddings;
  }
}

// /**
//  * Computes the lagged values for each element in a grid matrix based on a specified lag number and Moore neighborhood.
//  * For each element in the matrix, the function calculates the values of its neighbors at a specified lag distance
//  * in each of the 8 directions of the Moore neighborhood. If a neighbor is out of bounds, it is assigned a NaN value.
//  *
//  * Parameters:
//  *   mat    - A 2D vector representing the grid data.
//  *   lagNum - The number of steps to lag when considering the neighbors in the Moore neighborhood.
//  *
//  * Returns:
//  *   A 2D vector containing the lagged values for each element in the grid, arranged by the specified lag number.
//  *   If a neighbor is out of bounds, it is filled with NaN.
//  *
//  * Note:
//  *   The return value for each element is the lagged value of the neighbors, not the index of the neighbor.
//  */
// std::vector<std::vector<double>> CppLaggedVar4Grid(
//     const std::vector<std::vector<double>>& mat,
//     int lagNum
// ) {
//   int numRows = mat.size();
//   int numCols = mat[0].size();
//   int totalElements = numRows * numCols;
//
//   // Allocate space for each element, with 8 * lagNum neighbors, initially filled with NaN
//   std::vector<std::vector<double>> result(totalElements, std::vector<double>(8 * lagNum, std::numeric_limits<double>::quiet_NaN()));
//
//   // Directions for the Moore neighborhood (8 directions)
//   const int dx[] = {-1, -1, -1,  0, 0, 1, 1, 1};
//   const int dy[] = {-1,  0,  1, -1, 1,-1, 0, 1};
//
//   // Iterate over each element in the grid
//   for (int r = 0; r < numRows; ++r) {
//     for (int c = 0; c < numCols; ++c) {
//       int elementIndex = r * numCols + c;
//
//       // For each lag (from 1 to lagNum)
//       for (int lag = 1; lag <= lagNum; ++lag) {
//         // For each direction (8 directions)
//         for (int dir = 0; dir < 8; ++dir) {
//           // Calculate the correct index for storing the lagged values
//           int resultIndex = (lag - 1) * 8 + dir;
//           int dr = dx[dir] * lag; // row displacement for this direction
//           int dc = dy[dir] * lag; // column displacement for this direction
//           int rNeighbor = r + dr; // row of the neighbor at lag distance
//           int cNeighbor = c + dc; // column of the neighbor at lag distance
//
//           // Check if the neighbor position is within bounds of the grid
//           if (rNeighbor >= 0 && rNeighbor < numRows && cNeighbor >= 0 && cNeighbor < numCols) {
//             // Assign the value of the neighbor to the result matrix
//             result[elementIndex][resultIndex] = mat[rNeighbor][cNeighbor];
//           }
//           // If out of bounds, the value remains NaN (no need to assign anything)
//         }
//       }
//     }
//   }
//
//   // Return the result matrix containing lagged values for each element in the grid
//   return result;
// }

// /**
//  * Generates grid embeddings by calculating lagged variables for each element in a grid matrix,
//  * and stores the results in a matrix where each row represents an element and each column represents
//  * a different lagged value or the original element.
//  *
//  * Parameters:
//  *   mat  - A 2D vector representing the grid data.
//  *   E    - The number of embedding dimensions (columns in the resulting matrix).
//  *   tau  - The spatial lag step for constructing lagged state-space vectors.
//  *
//  * Returns:
//  *   A 2D vector (matrix) where each row contains the averaged lagged variables for
//  *   each embedding dimension (column).
//  */
// std::vector<std::vector<double>> GenGridEmbeddings(
//     const std::vector<std::vector<double>>& mat,
//     int E,
//     int tau
// ) {
//   int numRows = mat.size();
//   int numCols = mat[0].size();
//   int total_elements = numRows * numCols;
//
//   // Initialize the result matrix with total_elements rows and E columns
//   // Initially fill with NaN values
//   std::vector<std::vector<double>> result(total_elements, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));
//
//   if (tau == 0) {
//     // Flatten the matrix (mat) into the first column of the result matrix
//     int row = 0;
//     for (const auto& subset : mat) {
//       for (double value : subset) {
//         result[row][0] = value;
//         ++row;
//       }
//     }
//
//     // Fill the remaining columns (2 to E) with the averaged lagged variables
//     for (int lagNum = 1; lagNum < E; ++lagNum) {
//       // Calculate the lagged variables for the current lagNum
//       std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);
//
//       // Fill the current column (lagNum) with the averaged lagged variables
//       row = 0;
//       for (const auto& subset : lagged_vars) {
//         double sum = 0.0;
//         int count = 0;
//         for (double val : subset) {
//           if (!std::isnan(val)) {
//             sum += val;
//             ++count;
//           }
//         }
//
//         if (count > 0) {
//           result[row][lagNum] = sum / count; // Average the valid values
//         }
//         ++row;
//       }
//     }
//   } else {
//     // When tau != 0, start filling the result matrix from the tau-th column
//     int row = 0;
//     for (int lagNum = tau; lagNum < E + tau; ++lagNum) {
//       // Calculate the lagged variables for the current lagNum
//       std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);
//
//       // Fill the current column (lagNum - tau) with the averaged lagged variables
//       row = 0;
//       for (const auto& subset : lagged_vars) {
//         double sum = 0.0;
//         int count = 0;
//         for (double val : subset) {
//           if (!std::isnan(val)) {
//             sum += val;
//             ++count;
//           }
//         }
//
//         if (count > 0) {
//           result[row][lagNum - tau] = sum / count; // Average the valid values
//         }
//         ++row;
//       }
//     }
//   }
//
//   // Return the result matrix with grid embeddings
//   return result;
// }

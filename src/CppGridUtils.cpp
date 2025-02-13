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
int LocateGridIndices(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
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
  if (mat.empty() || mat[0].empty() || lagNum < 1) {
    return {};
  }

  const int rows = mat.size();
  const int cols = mat[0].size();
  const int numCells = rows * cols;
  const int numNeighbors = 8 * lagNum;

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
*   A 2D vector (matrix) where each row contains the original value (if includeself is true)
*   and the averaged lagged variables for each embedding dimension (column).
*
* If includeself is true, the first column will contain the original values from mat,
* and the subsequent columns will contain averaged lagged variables computed using the specified lag numbers.
* If includeself is false, the matrix will only contain the averaged lagged variables.
*/
std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    int tau
) {
  int numRows = mat.size();
  int numCols = mat[0].size();
  int total_elements = numRows * numCols;

  // Initialize the result matrix with total_elements rows and E columns
  // Initially fill with NaN values
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
    // When tau != 0, start filling the result matrix from the tau-th column
    int row = 0;
    for (int lagNum = tau; lagNum < E + tau; ++lagNum) {
      // Calculate the lagged variables for the current lagNum
      std::vector<std::vector<double>> lagged_vars = CppLaggedVar4Grid(mat, lagNum);

      // Fill the current column (lagNum - tau) with the averaged lagged variables
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
          result[row][lagNum - tau] = sum / count; // Average the valid values
        }
        ++row;
      }
    }
  }

  // Return the result matrix with grid embeddings
  return result;
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
//  *   A 2D vector (matrix) where each row contains the original value (if includeself is true)
//  *   and the averaged lagged variables for each embedding dimension (column).
//  *
//  * If includeself is true, the first column will contain the original values from mat,
//  * and the subsequent columns will contain averaged lagged variables computed using the specified lag numbers.
//  * If includeself is false, the matrix will only contain the averaged lagged variables.
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

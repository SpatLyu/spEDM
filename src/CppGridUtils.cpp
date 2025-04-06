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
    size_t k){
  const size_t rows = mat.size();
  const size_t cols = mat[0].size();
  std::vector<double> result(rows * cols, std::numeric_limits<double>::quiet_NaN());

  std::vector<double> me(rows * cols, std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      me.push_back(mat[i][j]);
    }
  }
  // The median of series me
  double s_me = CppMedian(me, true);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      const double center = mat[i][j];
      if (std::isnan(center)) continue;

      // Collect valid neighbors with absolute differences
      std::vector<std::pair<double, double>> neighbors; // <abs_diff, value>

      int max_lag = std::max({
        static_cast<int>(i),
        static_cast<int>(rows-1 - i),
        static_cast<int>(j),
        static_cast<int>(cols-1 - j)
      });

      // Expand lag until enough neighbors or reach max possible lag
      for (int lag = 1; lag <= max_lag; ++lag) {
        // Generate queen's case offsets for current lag
        std::vector<std::pair<int, int>> offsets;
        for (int dx = -lag; dx <= lag; ++dx) {
          for (int dy = -lag; dy <= lag; ++dy) {
            if (std::max(std::abs(dx), std::abs(dy)) == lag) {
              int ni = i + dx;
              int nj = j + dy;
              if (ni >= 0 && ni < static_cast<int>(rows) && nj >= 0 && nj < static_cast<int>(cols)) {
                double val = mat[ni][nj];
                if (!std::isnan(val)) {
                  neighbors.emplace_back(std::abs(val - center), val);
                }
              }
            }
          }
        }
        if (neighbors.size() >= k) break;
      }

      // Select top-k smallest differences
      if (neighbors.size() >= k) {
        // Sort by absolute difference then value
        std::sort(neighbors.begin(), neighbors.end(),
             [](const auto& a, const auto& b) {
               return a.first != b.first ? a.first < b.first : a.second < b.second;
             });

        // Symbolization
        std::vector<double> selected;
        for (size_t idx = 0; idx < k; ++idx){
          selected.push_back(neighbors[idx].second);
        }

        // The first indicator function
        std::vector<double> tau_s(selected.size());
        for(size_t i = 0; i < selected.size(); ++i){
          double tau = 0;
          if (selected[i] >= s_me){
            tau = 1;
          }
          tau_s[i] = tau;
        }

        // The second indicator function and the symbolization map
        std::vector<double> l_s(k);
        double taus = center >= s_me ? 1 : 0;
        for(size_t i = 0; i < k; ++i){
          double l = 0;
          if (tau_s[i] == taus){
            l = 1;
          }
          l_s[i] = l;
        }

        double fs = 0.0;
        for (size_t i = 0; i < l_s.size(); ++i) {
          fs += l_s[i];
        }

        result[i * cols + j] = fs;
      }
    }
  }

  return result;
}

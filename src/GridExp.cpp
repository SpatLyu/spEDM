#include <vector>
#include <cmath>
#include <algorithm>
#include "CppStats.h"
#include "CppGridUtils.h"
#include "Forecast4Grid.h"
#include "MultiViewEmbedding.h"
#include "GCCM4Grid.h"
#include "SCPCM4Grid.h"
#include "CrossMappingCardinality.h"
#include "SCT4Grid.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>

// [[Rcpp::export]]
int RcppLocateGridIndices(int curRow, int curCol,
                          int totalRow, int totalCol){
  int indices = LocateGridIndices(curRow,curCol,totalRow,totalCol);
  return indices + 1;
};

// [[Rcpp::export]]
Rcpp::NumericVector RcppRowColFromGrid(int cellNum, int totalCol){
  std::vector<int> result = RowColFromGrid(cellNum - 1, totalCol);
  for (int& val : result) {
    val += 1;
  }
  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
};

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppLaggedVar4Grid(const Rcpp::NumericMatrix& mat, int lagNum) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the CppLaggedVar4Grid function
  std::vector<std::vector<double>> laggedMat = CppLaggedVar4Grid(cppMat, lagNum);

  // Convert the result back to Rcpp::NumericMatrix
  int laggedRows = laggedMat.size();
  int laggedCols = laggedMat[0].size();
  Rcpp::NumericMatrix result(laggedRows, laggedCols);

  for (int r = 0; r < laggedRows; ++r) {
    for (int c = 0; c < laggedCols; ++c) {
      result(r, c) = laggedMat[r][c];
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGenGridEmbeddings(const Rcpp::NumericMatrix& mat,
                                          int E, int tau) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the GenGridEmbeddings function
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E, tau);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  int rows = embeddings.size();
  int cols = embeddings[0].size();
  Rcpp::NumericMatrix result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = embeddings[i][j];
    }
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector RcppGenGridSymbolization(const Rcpp::NumericMatrix& mat,
                                             int k) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the GenGridSymbolization function
  std::vector<double> symbolmap = GenGridSymbolization(cppMat, static_cast<size_t>(std::abs(k)));

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(symbolmap);
}

// [[Rcpp::export]]
Rcpp::IntegerVector RcppDivideGrid(const Rcpp::NumericMatrix& mat,
                                   int b, int shape = 3) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the CppDivideGrid function
  std::vector<int> blocks = CppDivideGrid(cppMat, b, shape);

  // Convert the result back to Rcpp::IntegerVector
  return Rcpp::wrap(blocks);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& mat,
                                     const Rcpp::IntegerMatrix& lib,
                                     const Rcpp::IntegerMatrix& pred,
                                     const Rcpp::IntegerVector& E,
                                     const Rcpp::IntegerVector& b,
                                     int tau,
                                     int threads) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> pred_indices(numRows * numCols, false);
  std::vector<bool> lib_indices(numRows * numCols, false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_indices[lib(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      if (!std::isnan(cppMat[currow-1][curcol-1])){
        lib_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_indices[pred(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      if (!std::isnan(cppMat[currow-1][curcol-1])){
        pred_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  std::vector<std::vector<double>> res_std = Simplex4Grid(
    cppMat,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "rho", "mae", "rmse");
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSMap4Grid(const Rcpp::NumericMatrix& mat,
                                  const Rcpp::IntegerMatrix& lib,
                                  const Rcpp::IntegerMatrix& pred,
                                  const Rcpp::NumericVector& theta,
                                  int E,
                                  int tau,
                                  int b,
                                  int threads) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> pred_indices(numRows * numCols, false);
  std::vector<bool> lib_indices(numRows * numCols, false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_indices[lib(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      if (!std::isnan(cppMat[currow-1][curcol-1])){
        lib_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_indices[pred(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      if (!std::isnan(cppMat[currow-1][curcol-1])){
        pred_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  // Convert Rcpp::NumericMatrix to std::vector<double>
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  std::vector<std::vector<double>> res_std = SMap4Grid(
    cppMat,
    lib_indices,
    pred_indices,
    theta_std,
    E,
    tau,
    b,
    threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("theta", "rho", "mae", "rmse");
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppMultiView4Grid(const Rcpp::NumericMatrix& xMatrix,
                                       const Rcpp::NumericMatrix& yMatrix,
                                       const Rcpp::IntegerMatrix& lib,
                                       const Rcpp::IntegerMatrix& pred,
                                       int E,
                                       int tau,
                                       int b,
                                       int top,
                                       int nvar,
                                       int threads){
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }
  std::vector<double> target = GridMat2Vec(yMatrix_cpp);

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> pred_indices(numRows * numCols, false);
  std::vector<bool> lib_indices(numRows * numCols, false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_indices[lib(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
        lib_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_indices[pred(i,0)-1] = true;
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
        pred_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
      }
    }
  }

  int num_row = xMatrix.nrow();
  int num_var = xMatrix.ncol();

  //  if top <= 0, we choose to apply the heuristic of k (sqrt(m))
  int k;
  if (top <= 0){
    double m = CppCombine(num_var*E,nvar) - CppCombine(num_var*(E - 1),nvar);
    k = std::floor(std::sqrt(m));
  } else {
    k = top;
  }

  // Combine all the lags in the embeddings
  std::vector<std::vector<double>> vec_std(num_row,std::vector<double>(E*num_var,std::numeric_limits<double>::quiet_NaN()));
  for (int n = 0; n < num_var; ++n) {
    // Initialize a std::vector to store the column values
    std::vector<double> univec(num_row);

    // Copy the nth column from the matrix to the vector
    for (int i = 0; i < num_row; ++i) {
      univec[i] = xMatrix(i, n);  // Access element at (i, n)
    }

    std::vector<std::vector<double>> unimat = GridVec2Mat(univec,numRows);

    // Generate the embedding:
    std::vector<std::vector<double>> vectors = GenGridEmbeddings(unimat,E,tau);

    for (size_t row = 0; row < vectors.size(); ++row) {  // Loop through each row
      for (size_t col = 0; col < vectors[0].size(); ++col) {  // Loop through each column
        vec_std[row][n * E + col] = vectors[row][col];  // Copy elements
      }
    }
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns; // To store indices of valid columns

  // Iterate over each column to check if it contains any non-NaN values
  for (size_t col = 0; col < vec_std[0].size(); ++col) {
    bool isAllNaN = true;
    for (size_t row = 0; row < vec_std.size(); ++row) {
      if (!std::isnan(vec_std[row][col])) {
        isAllNaN = false;
        break;
      }
    }
    if (!isAllNaN) {
      validColumns.push_back(col); // Store the index of valid columns
    }
  }

  if (validColumns.size() != vec_std[0].size()) {
    std::vector<std::vector<double>> filteredEmbeddings;
    for (size_t row = 0; row < vec_std.size(); ++row) {
      std::vector<double> filteredRow;
      for (size_t col : validColumns) {
        filteredRow.push_back(vec_std[row][col]);
      }
      filteredEmbeddings.push_back(filteredRow);
    }
    vec_std = filteredEmbeddings;
  }

  std::vector<double> res = MultiViewEmbedding(
    vec_std,
    target,
    lib_indices,
    pred_indices,
    b,
    k,
    threads);

  // Initialize a NumericMatrix with the given dimensions
  Rcpp::NumericMatrix resmat(numRows, numCols);

  // Fill the matrix with elements from res (assuming row-major order)
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      resmat(i, j) = res[i * numCols + j];  // Access element in row-major order
    }
  }

  return resmat;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    bool progressbar) {
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert libsizes to a fundamental C++ data type
  int libsizes_dim = libsizes.ncol();
  std::vector<int> libsizes_cpp1;
  std::vector<std::vector<int>> libsizes_cpp2(2);
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_cpp1.push_back(libsizes(i, 0));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_cpp2[0].push_back(libsizes(i, 0));
      libsizes_cpp2[1].push_back(libsizes(i, 1));
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_cpp1;
  std::vector<std::pair<int, int>> lib_cpp2(lib.nrow());
  if (libsizes_dim == 1){
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        lib_cpp1.push_back(lib(i, 0));
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        lib_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        lib_cpp2[i] = std::make_pair(lib(i, 0), lib(i, 1));
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<int> pred_cpp1;
  std::vector<std::pair<int, int>> pred_cpp2(pred.nrow());
  if (libsizes_dim == 1){
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        pred_cpp1.push_back(pred(i, 0));
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        pred_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        pred_cpp2[i] = std::make_pair(pred(i, 0), pred(i, 1));
      }
    }
  }

  // Call the C++ function GCCM4Grid or GCCM4GridOneDim
  std::vector<std::vector<double>> result;
  if (libsizes_dim == 1){
    result = GCCM4GridOneDim(
      xMatrix_cpp,
      yMatrix_cpp,
      libsizes_cpp1,
      lib_cpp1,
      pred_cpp1,
      E,
      tau,
      b,
      simplex,
      theta,
      threads,
      parallel_level,
      progressbar
    );
  } else{
    result = GCCM4Grid(
      xMatrix_cpp,
      yMatrix_cpp,
      libsizes_cpp2,
      lib_cpp2,
      pred_cpp2,
      E,
      tau,
      b,
      simplex,
      theta,
      threads,
      parallel_level,
      progressbar
    );
  }

  Rcpp::NumericMatrix resultMatrix(result.size(), 5);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("libsizes",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_upper","x_xmap_y_lower");
  return resultMatrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSCPCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::NumericMatrix& zMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    const Rcpp::IntegerVector& b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    bool cumulate,
    bool progressbar) {
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> zMatrix_cpp(zMatrix.ncol());
  for (int i = 0; i < zMatrix.ncol(); ++i) {
    Rcpp::NumericVector covvar = zMatrix.column(i);
    zMatrix_cpp[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert libsizes to a fundamental C++ data type
  int libsizes_dim = libsizes.ncol();
  std::vector<int> libsizes_cpp1;
  std::vector<std::vector<int>> libsizes_cpp2(2);
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_cpp1.push_back(libsizes(i, 0));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_cpp2[0].push_back(libsizes(i, 0));
      libsizes_cpp2[1].push_back(libsizes(i, 1));
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_cpp1;
  std::vector<std::pair<int, int>> lib_cpp2(lib.nrow());
  if (libsizes_dim == 1){
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        lib_cpp1.push_back(lib(i, 0));
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        lib_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        lib_cpp2[i] = std::make_pair(lib(i, 0), lib(i, 1));
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<int> pred_cpp1;
  std::vector<std::pair<int, int>> pred_cpp2(pred.nrow());
  if (libsizes_dim == 1){
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        pred_cpp1.push_back(pred(i, 0));
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        if (!std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        pred_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        pred_cpp2[i] = std::make_pair(pred(i, 0), pred(i, 1));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_cpp = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_cpp = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_cpp = Rcpp::as<std::vector<int>>(b);

  // Call the C++ function SCPCM4Grid or SCPCM4GridOneDim
  std::vector<std::vector<double>> result;
  if (libsizes_dim == 1){
    result = SCPCM4GridOneDim(
      xMatrix_cpp,
      yMatrix_cpp,
      zMatrix_cpp,
      libsizes_cpp1,
      lib_cpp1,
      pred_cpp1,
      E_cpp,
      tau_cpp,
      b_cpp,
      simplex,
      theta,
      threads,
      parallel_level,
      cumulate,
      progressbar
    );
  } else{
    result = SCPCM4Grid(
      xMatrix_cpp,
      yMatrix_cpp,
      zMatrix_cpp,
      libsizes_cpp2,
      lib_cpp2,
      pred_cpp2,
      E_cpp,
      tau_cpp,
      b_cpp,
      simplex,
      theta,
      threads,
      parallel_level,
      cumulate,
      progressbar
    );
  }

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 9);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
    resultMatrix(i, 5) = result[i][5];
    resultMatrix(i, 6) = result[i][6];
    resultMatrix(i, 7) = result[i][7];
    resultMatrix(i, 8) = result[i][8];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create(
    "libsizes","T_mean","D_mean",
    "T_sig","T_upper","T_lower",
    "D_sig","D_upper","D_lower");
  return resultMatrix;
}

// Wrapper function to perform GCMC for spatial grid data
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCMC4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    const Rcpp::IntegerVector& b,
    const Rcpp::IntegerVector& max_r,
    int threads,
    bool progressbar){
  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  double validCellNum = 0;
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      if (!std::isnan(yMatrix(i, j))){
        validCellNum += 1;
      }
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);
  std::vector<int> maxr_std = Rcpp::as<std::vector<int>>(max_r);

  // Remove values in b_std that are greater than validCellNum or less than or equal to 3
  b_std.erase(std::remove_if(b_std.begin(), b_std.end(),
                             [validCellNum](int x) { return x > validCellNum || x <= 3; }),
                             b_std.end());

  if (b_std.empty()) {
    Rcpp::stop("k cannot be less than or equal to 3 or greater than the number of non-NA values.");
  }

  // Remove duplicates for b_std
  std::sort(b_std.begin(), b_std.end());
  b_std.erase(std::unique(b_std.begin(), b_std.end()), b_std.end());

  // Convert Rcpp IntegerMatrix to std::vector<int>
  int n_libcol = lib.ncol();
  int n_predcol = lib.ncol();
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  std::vector<int> lib_std;
  if (n_libcol == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std.push_back(lib(i,0) - 1);
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int rowLibIndice = lib(i,0);
      int colLibIndice = lib(i,1);
      if (!std::isnan(yMatrix_cpp[rowLibIndice-1][colLibIndice-1])){
        lib_std.push_back(LocateGridIndices(rowLibIndice, colLibIndice, numRows, numCols));
      }
    }
  }

  std::vector<int> pred_std;
  if (n_predcol == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_std.push_back(pred(i,0) - 1);
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int rowPredIndice = pred(i,0);
      int colPredIndice = pred(i,1);
      if (!std::isnan(yMatrix_cpp[rowPredIndice-1][colPredIndice-1])){
        pred_std.push_back(LocateGridIndices(rowPredIndice, colPredIndice, numRows, numCols));
      }
    }
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenGridEmbeddings(xMatrix_cpp, E[0], tau_std[0]);
  std::vector<std::vector<double>> e2 = GenGridEmbeddings(yMatrix_cpp, E[1], tau_std[1]);

  // Perform GCMC for spatial grid data
  std::vector<std::vector<double>> cs1 = CrossMappingCardinality(e1,e2,lib_std,pred_std,b_std,maxr_std,threads,progressbar);

  Rcpp::NumericMatrix resultMatrix(b_std.size(), 5);
  for (size_t i = 0; i < b_std.size(); ++i) {
    for (size_t j = 0; j < cs1[0].size(); ++j){
      resultMatrix(i, j) = cs1[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("neighbors",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_upper","x_xmap_y_lower");
  return resultMatrix;
}

// Wrapper function to perform SCT for spatial grid data
// [[Rcpp::export]]
Rcpp::NumericVector RcppSCT4Grid(const Rcpp::NumericMatrix& x,
                                 const Rcpp::NumericMatrix& y,
                                 const Rcpp::IntegerMatrix& block,
                                 int k,
                                 int threads,
                                 int boot = 399,
                                 double base = 2,
                                 unsigned int seed = 42,
                                 bool symbolize = true,
                                 bool normalize = false,
                                 bool progressbar = true){
  int numRows = y.nrow();
  int numCols = y.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xmat(x.nrow(), std::vector<double>(x.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      xmat[i][j] = x(i, j);
    }
  }
  std::vector<std::vector<double>> ymat(y.nrow(), std::vector<double>(y.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      ymat[i][j] = y(i, j);
    }
  }

  // Convert block to a fundamental C++ data type
  int b_dim = block.ncol();
  std::vector<int> b_std;
  if (b_dim == 1){
    for (int i = 0; i < block.nrow(); ++i) {
      b_std.push_back(block(i, 0));
    }
  } else {
    for (int i = 0; i < block.nrow(); ++i) {
        // Convert to 0-based index
        int currow = block(i,0);
        int curcol = block(i,1);
        b_std.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
  }

  // Perform SCT for spatial grid data
  std::vector<double> sc = SCT4Grid(
    xmat,
    ymat,
    b_std,
    k,
    threads,
    boot,
    base,
    seed,
    symbolize,
    normalize,
    progressbar
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "significance for x → y causality",
    "statistic for y → x causality",
    "significance for y → x causality"
  );

  return sc_res;
}

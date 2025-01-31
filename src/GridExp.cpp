#include <vector>
#include "CppGridUtils.h"
#include "Forecast4Grid.h"
#include "GCCM4Grid.h"
#include "SCPCM4Grid.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppLaggedVar4Grid(Rcpp::NumericMatrix mat, int lagNum) {
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
Rcpp::NumericMatrix RcppGenGridEmbeddings(Rcpp::NumericMatrix mat, int E, bool includeself) {
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
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E, includeself);

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
int RcppLocateGridIndices(int curRow, int curCol,
                          int totalRow, int totalCol){
  int indices = LocateGridIndices(curRow,curCol,totalRow,totalCol);
  return indices + 1;
};

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& mat,
                                     const Rcpp::IntegerMatrix& lib,
                                     const Rcpp::IntegerMatrix& pred,
                                     const Rcpp::IntegerVector& E,
                                     int b,
                                     int threads,
                                     bool includeself) {
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

  for (int i = 0; i < lib.nrow(); ++i) {
    // Convert to 0-based index
    currow = lib(i,0);
    curcol = lib(i,1);
    if (!std::isnan(cppMat[currow-1][curcol-1])){
      lib_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
    }
  }

  for (int i = 0; i < pred.nrow(); ++i) {
    // Convert to 0-based index
    currow = pred(i,0);
    curcol = pred(i,1);
    if (!std::isnan(cppMat[currow-1][curcol-1])){
      pred_indices[LocateGridIndices(currow, curcol, numRows, numCols)] = true;
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  std::vector<std::vector<double>> res_std = Simplex4Grid(
    cppMat,
    lib_indices,
    pred_indices,
    E_std,
    b,
    threads,
    includeself);

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
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "rho", "mae", "rmse");
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGCCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerVector& lib_sizes,
    const Rcpp::IntegerMatrix& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    bool includeself,
    bool progressbar) {

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

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_sizes_cpp(lib_sizes.size());
  for (int i = 0; i < lib_sizes.size(); ++i) {
    lib_sizes_cpp[i] = lib_sizes[i];
  }

  // Convert Rcpp IntegerMatrix to std::vector<std::pair<int, int>>
  std::vector<std::pair<int, int>> pred_cpp(pred.nrow());
  for (int i = 0; i < pred.nrow(); ++i) {
    pred_cpp[i] = std::make_pair(pred(i, 0), pred(i, 1));
  }

  // Call the C++ function GCCM4Grid
  std::vector<std::vector<double>> result = GCCM4Grid(
    xMatrix_cpp,
    yMatrix_cpp,
    lib_sizes_cpp,
    pred_cpp,
    E,
    tau,
    b,
    simplex,
    theta,
    threads,
    includeself,
    progressbar
  );

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
    const Rcpp::NumericMatrix& zMatrixs,
    const Rcpp::IntegerVector& lib_sizes,
    const Rcpp::IntegerMatrix& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    bool includeself,
    bool progressbar) {

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

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> lib_sizes_cpp(lib_sizes.size());
  for (int i = 0; i < lib_sizes.size(); ++i) {
    lib_sizes_cpp[i] = lib_sizes[i];
  }

  // Convert Rcpp IntegerMatrix to std::vector<std::pair<int, int>>
  std::vector<std::pair<int, int>> pred_cpp(pred.nrow());
  for (int i = 0; i < pred.nrow(); ++i) {
    pred_cpp[i] = std::make_pair(pred(i, 0), pred(i, 1));
  }

  // Call the C++ function GCCM4Grid
  std::vector<std::vector<double>> result = GCCM4Grid(
    xMatrix_cpp,
    yMatrix_cpp,
    lib_sizes_cpp,
    pred_cpp,
    E,
    tau,
    b,
    simplex,
    theta,
    threads,
    includeself,
    progressbar
  );

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

// // [[Rcpp::export]]
// Rcpp::List RcppGenGridEmbeddings2(Rcpp::NumericMatrix mat, int E) {
//   // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
//   int numRows = mat.nrow();
//   int numCols = mat.ncol();
//   std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));
//
//   for (int r = 0; r < numRows; ++r) {
//     for (int c = 0; c < numCols; ++c) {
//       cppMat[r][c] = mat(r, c);
//     }
//   }
//
//   // Call the GenGridEmbeddings function
//   std::vector<std::vector<std::vector<double>>> embeddings = GenGridEmbeddings2(cppMat, E);
//
//   // Convert the result back to an Rcpp::List of Rcpp::NumericMatrix
//   Rcpp::List result(E + 1);
//
//   for (int i = 0; i <= E; ++i) {
//     int embeddingRows = embeddings[i].size();
//     int embeddingCols = embeddings[i][0].size();
//     Rcpp::NumericMatrix embeddingMat(embeddingRows, embeddingCols);
//
//     for (int r = 0; r < embeddingRows; ++r) {
//       for (int c = 0; c < embeddingCols; ++c) {
//         embeddingMat(r, c) = embeddings[i][r][c];
//       }
//     }
//
//     result[i] = embeddingMat;
//   }
//
//   return result;
// }

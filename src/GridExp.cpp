#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "CppGridUtils.h"
#include "Forecast4Grid.h"
#include "MultiViewEmbedding.h"
#include "GCCM4Grid.h"
#include "SCPCM4Grid.h"
#include "CrossMappingCardinality.h"
#include "FalseNearestNeighbors.h"
#include "SLM4Grid.h"
#include "SGC4Grid.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Wrapper function to convert row and column indices to 1-based linear indices
// [[Rcpp::export(rng = false)]]
int RcppLocateGridIndices(int curRow, int curCol,
                          int totalRow, int totalCol){
  int indices = LocateGridIndices(curRow,curCol,totalRow,totalCol);
  return indices + 1;
};

// Wrapper function to convert a row-major linear index to row and column indices
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppRowColFromGrid(int cellNum, int totalCol){
  std::vector<int> result = RowColFromGrid(cellNum - 1, totalCol);
  for (int& val : result) {
    val += 1;
  }
  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
};

// Wrapper function to calculate lagged values for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppLaggedVal4Grid(const Rcpp::NumericMatrix& mat, int lagNum = 1) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the CppLaggedVal4Grid function
  std::vector<std::vector<double>> laggedMat = CppLaggedVal4Grid(cppMat, lagNum);

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

// Wrapper function to generate embeddings for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenGridEmbeddings(const Rcpp::NumericMatrix& mat,
                                          int E = 3, int tau = 1, int style = 1) {
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
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E, tau, style);

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

// Wrapper function to generate neighbors for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGenGridNeighbors(const Rcpp::NumericMatrix& mat,
                                const Rcpp::IntegerMatrix& lib,
                                int k = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_std(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std[i] = lib(i, 0) - 1;
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std[i] = LocateGridIndices(lib(i, 0), lib(i, 1), numRows, numCols);
    }
  }

  // Generate neighbors
  std::vector<std::vector<int>> neighbors = GenGridNeighbors(
    cppMat, lib_std, static_cast<size_t>(std::abs(k))
  );

  // Convert neighbors to Rcpp::List with 1-based indexing
  int n = neighbors.size();
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int> neighbor_i = neighbors[i];
    for (auto& idx : neighbor_i) {
      idx += 1; // convert to 1-based index
    }
    result[i] = Rcpp::IntegerVector(neighbor_i.begin(), neighbor_i.end());
  }

  return result;
}

// Wrapper function to implement a symbolic transformation of a univariate spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenGridSymbolization(const Rcpp::NumericMatrix& mat,
                                             const Rcpp::IntegerMatrix& lib,
                                             const Rcpp::IntegerMatrix& pred,
                                             int k = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      lib_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std[i] = std::make_pair(lib(i, 0) - 1, lib(i, 1) - 1);
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      pred_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_std[i] = std::make_pair(pred(i, 0) - 1, pred(i, 1) - 1);
    }
  }

  // Call the GenGridSymbolization function
  std::vector<double> symbolmap = GenGridSymbolization(
    cppMat, lib_std, pred_std, static_cast<size_t>(std::abs(k))
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(symbolmap);
}

// Wrapper function to partition spatial units in spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppDivideGrid(const Rcpp::NumericMatrix& mat,
                                   int b = 9, int shape = 3) {
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

// Wrapper function to perform univariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSLMUni4Grid(
    const Rcpp::NumericMatrix& mat,
    int k = 4,
    int step = 20,
    double alpha = 0.77,
    double escape_threshold = 1e10
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<double>> result = SLMUni4Grid(cppMat, k, step, alpha, escape_threshold);

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result.size());
  int n_cols = result[0].size();
  Rcpp::NumericMatrix out(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out(i, j) = result[i][j];
    }
  }

  return out;
}

// Wrapper function to perform bivariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMBi4Grid(
    const Rcpp::NumericMatrix& mat1,
    const Rcpp::NumericMatrix& mat2,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double beta_xy = 0.05,
    double beta_yx = 0.4,
    double escape_threshold = 1e10
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat1.nrow();
  int numCols = mat1.ncol();
  std::vector<std::vector<double>> cppMat1(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat2(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat1[r][c] = mat1(r, c);
      cppMat2[r][c] = mat2(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMBi4Grid(
    cppMat1, cppMat2, k, step, alpha_x, alpha_y, beta_xy, beta_yx, escape_threshold
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = result[0][1].size();
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y
  );

  return out;
}

// Wrapper function to perform trivariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMTri4Grid(
    const Rcpp::NumericMatrix& mat1,
    const Rcpp::NumericMatrix& mat2,
    const Rcpp::NumericMatrix& mat3,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double alpha_z = 0.55,
    double beta_xy = 0.05,
    double beta_xz = 0.05,
    double beta_yx = 0.4,
    double beta_yz = 0.4,
    double beta_zx = 0.65,
    double beta_zy = 0.65,
    double escape_threshold = 1e10
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat1.nrow();
  int numCols = mat1.ncol();
  std::vector<std::vector<double>> cppMat1(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat2(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat3(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat1[r][c] = mat1(r, c);
      cppMat2[r][c] = mat2(r, c);
      cppMat3[r][c] = mat3(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMTri4Grid(
    cppMat1, cppMat2, cppMat3,
    k, step, alpha_x, alpha_y, alpha_z,
    beta_xy, beta_xz, beta_yx, beta_yz, beta_zx, beta_zy,
    escape_threshold
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = result[0][1].size();
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);
  Rcpp::NumericMatrix out_z(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
      out_z(i, j) = result[2][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y,
    Rcpp::Named("z") = out_z
  );

  return out;
}

// Wrapper function to perform FNN for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppFNN4Grid(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::NumericVector& rt,
    const Rcpp::NumericVector& eps,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    int tau = 1,
    int threads = 8,
    int parallel_level = 0){
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  double validCellNum = 0;
  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
      if (!std::isnan(mat(r, c))){
        validCellNum += 1;
      }
    }
  }

  // Convert Rcpp NumericVector to std::vector<double>
  std::vector<double> rt_std = Rcpp::as<std::vector<double>>(rt);
  std::vector<double> eps_std = Rcpp::as<std::vector<double>>(eps);

  // Convert Rcpp IntegerMatrix to std::vector<int>
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();

  std::vector<int> lib_std;
  if (n_libcol == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(cppMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_std.push_back(lib(i,0) - 1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int rowLibIndice = lib(i,0);
      int colLibIndice = lib(i,1);
      if (!std::isnan(cppMat[rowLibIndice-1][colLibIndice-1])){
        lib_std.push_back(LocateGridIndices(rowLibIndice, colLibIndice, numRows, numCols));
      }
    }
  }

  std::vector<int> pred_std;
  if (n_predcol == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(cppMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_std.push_back(pred(i,0) - 1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int rowPredIndice = pred(i,0);
      int colPredIndice = pred(i,1);
      if (!std::isnan(cppMat[rowPredIndice-1][colPredIndice-1])){
        pred_std.push_back(LocateGridIndices(rowPredIndice, colPredIndice, numRows, numCols));
      }
    }
  }

  // Generate embeddings
  std::vector<double> E_std = Rcpp::as<std::vector<double>>(E);
  int max_E = CppMax(E_std, true);
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, max_E, tau);

  // Perform FNN for spatial grid data
  std::vector<double> fnn = CppFNN(embeddings,lib_std,pred_std,rt_std,eps_std,true,threads,parallel_level);

  // Convert the result back to Rcpp::NumericVector and set names as "E:1", "E:2", ..., "E:n"
  Rcpp::NumericVector result = Rcpp::wrap(fnn);
  Rcpp::CharacterVector resnames(result.size());
  for (int i = 0; i < result.size(); ++i) {
    resnames[i] = "E:" + std::to_string(i + 1);
  }
  result.names() = resnames;

  return result;
}

// Wrapper function to perform simplex forecasting for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& source,
                                     const Rcpp::NumericMatrix& target,
                                     const Rcpp::IntegerMatrix& lib,
                                     const Rcpp::IntegerMatrix& pred,
                                     const Rcpp::IntegerVector& E,
                                     const Rcpp::IntegerVector& b,
                                     int tau = 1,
                                     int threads = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> pred_indices;
  std::vector<int> lib_indices;

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  std::vector<std::vector<double>> res_std = Simplex4Grid(
    sourceMat,
    targetMat,
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

// Wrapper function to perform s-mapping for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSMap4Grid(const Rcpp::NumericMatrix& source,
                                  const Rcpp::NumericMatrix& target,
                                  const Rcpp::IntegerMatrix& lib,
                                  const Rcpp::IntegerMatrix& pred,
                                  const Rcpp::NumericVector& theta,
                                  int E = 3,
                                  int tau = 1,
                                  int b = 5,
                                  int threads = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> pred_indices;
  std::vector<int> lib_indices;

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  // Convert Rcpp::NumericMatrix to std::vector<double>
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  std::vector<std::vector<double>> res_std = SMap4Grid(
    sourceMat,
    targetMat,
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

// Wrapper function to perform multiview forecasting for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppMultiView4Grid(const Rcpp::NumericMatrix& xMatrix,
                                       const Rcpp::NumericMatrix& yMatrix,
                                       const Rcpp::IntegerMatrix& lib,
                                       const Rcpp::IntegerMatrix& pred,
                                       int E = 3,
                                       int tau = 1,
                                       int b = 5,
                                       int top = 5,
                                       int nvar = 3,
                                       int threads = 8){
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(numRows, std::vector<double>(numCols));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert yMatrix to std::vector<double>
  std::vector<double> target(numRows * numCols);

  if (numCols == 1){
    for (int i = 0; i < numRows; ++i) {
      target[i] =  yMatrix(i,0);
    }
  } else {
    target = GridMat2Vec(yMatrix_cpp);
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> pred_indices;
  std::vector<int> lib_indices;

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      int currow = lib(i,0);
      int curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      int cellindice = LocateGridIndices(currow, curcol, numRows, numCols);
      if (!std::isnan(target[cellindice])){
        lib_indices.push_back(cellindice);
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      int currow = pred(i,0);
      int curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      int cellindice = LocateGridIndices(currow, curcol, numRows, numCols);
      if (!std::isnan(target[cellindice])){
        pred_indices.push_back(cellindice);
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

// Wrapper function to compute intersection cardinality for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppIC4Grid(const Rcpp::NumericMatrix& source,
                                const Rcpp::NumericMatrix& target,
                                const Rcpp::IntegerMatrix& lib,
                                const Rcpp::IntegerMatrix& pred,
                                const Rcpp::IntegerVector& E,
                                const Rcpp::IntegerVector& b,
                                int tau = 1,
                                int exclude = 0,
                                int threads = 8,
                                int parallel_level = 0) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<size_t> pred_indices;
  std::vector<size_t> lib_indices;

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(static_cast<size_t>(lib(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(static_cast<size_t>(pred(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  // Check the validity of the neignbor numbers
  std::vector<int> b_std;
  for (int i = 0; i < b.size(); ++i){
    if (b[i] > static_cast<int>(lib_indices.size())) {
      Rcpp::stop("Neighbor numbers count out of acceptable range at position %d (value: %d)", i + 1, b[i]);
    }
    b_std.push_back(b[i]);
  }

  std::vector<std::vector<double>> res_std = IC4Grid(
    sourceMat,
    targetMat,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    exclude,
    threads,
    parallel_level);

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
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "CausalScore", "Significance");
  return result;
}

// Wrapper function to perform GCCM for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGCCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    int E = 3,
    int tau = 1,
    int b = 5,
    bool simplex = true,
    double theta = 0,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = false) {
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
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
          lib_cpp1.push_back(lib(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          lib_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1])){
          lib_cpp2[i] = std::make_pair(lib(i, 0) - 1, lib(i, 1) - 1);
        }
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
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
          pred_cpp1.push_back(pred(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          pred_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1])){
          pred_cpp2[i] = std::make_pair(pred(i, 0) - 1, pred(i, 1) - 1);
        }
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

// Wrapper function to perform SCPCM for spatial grid data
// [[Rcpp::export(rng = false)]]
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
    bool simplex = true,
    double theta = 0,
    int threads = 8,
    int parallel_level = 0,
    bool cumulate = false,
    bool progressbar = false) {
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
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
          lib_cpp1.push_back(lib(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          lib_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1])){
          lib_cpp2[i] = std::make_pair(lib(i, 0) - 1, lib(i, 1) - 1);
        }
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
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
          pred_cpp1.push_back(pred(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          pred_cpp2[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1])){
          pred_cpp2[i] = std::make_pair(pred(i, 0) - 1, pred(i, 1) - 1);
        }
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
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGCMC4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int b,
    int r = 0,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = false){
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

  // Convert Rcpp IntegerMatrix to std::vector<int>
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert libsizes to a fundamental C++ data type
  int libsizes_dim = libsizes.ncol();
  std::vector<size_t> libsizes_std;
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_std.push_back(static_cast<size_t>(libsizes(i, 0)));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_std.push_back(static_cast<size_t>(LocateGridIndices(libsizes(i, 0),libsizes(i, 1),
                                                                   numRows, numCols)));
    }
  }

  std::vector<size_t> lib_std;
  if (n_libcol == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        lib_std.push_back(static_cast<size_t>(lib(i,0) - 1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int rowLibIndice = lib(i,0);
      int colLibIndice = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[rowLibIndice-1][colLibIndice-1]) &&
          !std::isnan(yMatrix_cpp[rowLibIndice-1][colLibIndice-1])){
        lib_std.push_back(static_cast<size_t>(LocateGridIndices(rowLibIndice, colLibIndice, numRows, numCols)));
      }
    }
  }

  std::vector<size_t> pred_std;
  if (n_predcol == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_std.push_back(static_cast<size_t>(pred(i, 0) - 1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int rowPredIndice = pred(i,0);
      int colPredIndice = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[rowPredIndice-1][colPredIndice-1]) &&
          !std::isnan(yMatrix_cpp[rowPredIndice-1][colPredIndice-1])){
        pred_std.push_back(static_cast<size_t>(LocateGridIndices(rowPredIndice, colPredIndice, numRows, numCols)));
      }
    }
  }

  if (b < 3 || b > validCellNum) {
    Rcpp::stop("k cannot be less than or equal to 3 or greater than the number of non-NA values.");
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenGridEmbeddings(xMatrix_cpp, E[0], tau_std[0]);
  std::vector<std::vector<double>> e2 = GenGridEmbeddings(yMatrix_cpp, E[1], tau_std[1]);

  // Perform GCMC for spatial grid data
  CMCRes res = CrossMappingCardinality(e1,e2,libsizes_std,lib_std,pred_std,
                                       static_cast<size_t>(b),static_cast<size_t>(r),
                                       threads,parallel_level,progressbar);

  // Convert mean_aucs to Rcpp::DataFrame
  std::vector<double> libs, aucs;
  for (const auto& cm : res.causal_strength) {
    libs.push_back(cm[0]);
    aucs.push_back(cm[1]);
  }

  Rcpp::DataFrame cs_df = Rcpp::DataFrame::create(
    Rcpp::Named("libsizes") = libs,
    Rcpp::Named("x_xmap_y_mean") = aucs
  );

  // Wrap causal_strength with names
  Rcpp::DataFrame xmap_df = Rcpp::DataFrame::create(
    Rcpp::Named("neighbors") = res.cross_mapping[0],
                                                Rcpp::Named("x_xmap_y_mean") = res.cross_mapping[1],
                                                                                                Rcpp::Named("x_xmap_y_sig") = res.cross_mapping[2],
                                                                                                                                               Rcpp::Named("x_xmap_y_upper") = res.cross_mapping[3],
                                                                                                                                                                                                Rcpp::Named("x_xmap_y_lower")  = res.cross_mapping[4]
  );

  return Rcpp::List::create(
    Rcpp::Named("xmap") = xmap_df,
    Rcpp::Named("cs") = cs_df
  );
}

// Wrapper function to perform SGC for spatial grid data without bootstrapped significance
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGCSingle4Grid(const Rcpp::NumericMatrix& x,
                                       const Rcpp::NumericMatrix& y,
                                       const Rcpp::IntegerMatrix& lib,
                                       const Rcpp::IntegerMatrix& pred,
                                       int k,
                                       double base = 2,
                                       bool symbolize = true,
                                       bool normalize = false){
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

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        lib_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      if (!std::isnan(xmat[lib(i, 0) - 1][lib(i, 1) - 1]) &&
          !std::isnan(ymat[lib(i, 0) - 1][lib(i, 1) - 1])){
        lib_std[i] = std::make_pair(lib(i, 0) - 1, lib(i, 1) - 1);
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        pred_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      if (!std::isnan(xmat[pred(i, 0) - 1][pred(i, 1) - 1]) &&
          !std::isnan(ymat[pred(i, 0) - 1][pred(i, 1) - 1])){
        pred_std[i] = std::make_pair(pred(i, 0) - 1, pred(i, 1) - 1);
      }
    }
  }

  // Perform SGC for spatial grid data
  std::vector<double> sc = SGCSingle4Grid(
    xmat,
    ymat,
    lib_std,
    pred_std,
    k,
    base,
    symbolize,
    normalize
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "statistic for y → x causality"
  );

  return sc_res;
}

// Wrapper function to perform SGC for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGC4Grid(const Rcpp::NumericMatrix& x,
                                 const Rcpp::NumericMatrix& y,
                                 const Rcpp::IntegerMatrix& lib,
                                 const Rcpp::IntegerMatrix& pred,
                                 const Rcpp::IntegerMatrix& block,
                                 int k,
                                 int threads = 8,
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

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        lib_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      if (!std::isnan(xmat[lib(i, 0) - 1][lib(i, 1) - 1]) &&
          !std::isnan(ymat[lib(i, 0) - 1][lib(i, 1) - 1])){
        lib_std[i] = std::make_pair(lib(i, 0) - 1, lib(i, 1) - 1);
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        pred_std[i] = std::make_pair(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      if (!std::isnan(xmat[pred(i, 0) - 1][pred(i, 1) - 1]) &&
          !std::isnan(ymat[pred(i, 0) - 1][pred(i, 1) - 1])){
        pred_std[i] = std::make_pair(pred(i, 0) - 1, pred(i, 1) - 1);
      }
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

  // Perform SGC for spatial grid data
  std::vector<double> sc = SGC4Grid(
    xmat,
    ymat,
    lib_std,
    pred_std,
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

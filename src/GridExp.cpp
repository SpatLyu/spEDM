#include <vector>
#include "CppGridUtils.h"
#include "GCCM4Grid.h"
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
Rcpp::NumericMatrix RcppGenGridEmbeddings(Rcpp::NumericMatrix mat, int E) {
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
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E);

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
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& mat,
                                     const Rcpp::IntegerMatrix& lib,
                                     const Rcpp::IntegerMatrix& pred,
                                     const Rcpp::IntegerVector& E,
                                     int b) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols); // Reserve space for efficiency

  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      vec_std.push_back(mat(i, j)); // Add element to the vector
    }
  }

  // Initialize lib_indices and pred_indices with all false
  std::vector<bool> pred_indices(numRows * numCols, false);
  std::vector<bool> lib_indices(numRows * numCols, false);

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  for (int i = 0; i < lib.nrow(); ++i) {
    lib_indices[LocateGridIndices(lib(i,0), lib(i,1), numRows, numCols)] = true; // Convert to 0-based index
  }
  for (int i = 0; i < pred.nrow(); ++i) {
    pred_indices[LocateGridIndices(lib(i,0), lib(i,1), numRows, numCols)] = true; // Convert to 0-based index
  }

  // Initialize the result matrix
  Rcpp::NumericMatrix result(E.size(), 4); // Columns: E, PearsonCor, MAE, RMSE

  // Loop over each embedding dimension E
  for (int i = 0; i < E.size(); ++i) {
    // Call the GenGridEmbeddings function
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E[i]);

    // Call SimplexBehavior to compute metrics
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, b);

    // Store results in the matrix
    result(i, 0) = E[i];               // Embedding dimension
    result(i, 1) = metrics[0];         // PearsonCor
    result(i, 2) = metrics[1];         // MAE
    result(i, 3) = metrics[2];         // RMSE
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

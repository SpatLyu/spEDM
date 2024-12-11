#include <vector>
#include "CppGridUtils.h"
#include <Rcpp.h>

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

  // Call the Cpp function
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
Rcpp::List RcppGenGridEmbeddings(Rcpp::NumericMatrix mat, int E) {
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
  std::vector<std::vector<std::vector<double>>> embeddings = GenGridEmbeddings(cppMat, E);

  // Convert the result back to an Rcpp::List of Rcpp::NumericMatrix
  Rcpp::List result(E + 1);

  for (int i = 0; i <= E; ++i) {
    int embeddingRows = embeddings[i].size();
    int embeddingCols = embeddings[i][0].size();
    Rcpp::NumericMatrix embeddingMat(embeddingRows, embeddingCols);

    for (int r = 0; r < embeddingRows; ++r) {
      for (int c = 0; c < embeddingCols; ++c) {
        embeddingMat(r, c) = embeddings[i][r][c];
      }
    }

    result[i] = embeddingMat;
  }

  return result;
}

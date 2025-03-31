#include <vector>
#include "Entropy.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double RcppEntropy(const Rcpp::NumericVector& vec,
                   int k = 3, double base = 10,
                   bool L1norm = false, bool NA_rm = false){
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec);
  return CppEntropy(v1,static_cast<size_t>(std::abs(k)),base,L1norm,NA_rm);
};

// [[Rcpp::export]]
double RcppJoinEntropy(const Rcpp::NumericMatrix& mat,
                       int k = 3, double base = 10,
                       bool L1norm = false, bool NA_rm = false){
  // Convert the Rcpp::NumericMatrix to a C++ vector of vectors (std::vector)
  size_t rownum = mat.nrow();
  size_t colnum = mat.ncol();
  std::vector<std::vector<double>> cppMat(rownum, std::vector<double>(colnum));

  // Fill cppMat with values from the R matrix
  for (size_t i = 0; i < rownum; ++i) {
    for (size_t j = 0; j < colnum; ++j) {
      cppMat[i][j] = mat(i, j);
    }
  }

  return CppJoinEntropy(cppMat,static_cast<size_t>(std::abs(k)),base,L1norm,NA_rm);
};

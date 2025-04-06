#include <vector>
#include "Entropy.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double RcppEntropy_Cont(const Rcpp::NumericVector& vec, int k = 3,
                   double base = 10, bool NA_rm = false){
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec);
  return CppEntropy_Cont(v1,static_cast<size_t>(std::abs(k)),base,NA_rm);
};

// [[Rcpp::export]]
double RcppJoinEntropy_Cont(const Rcpp::NumericVector& vec1,
                            const Rcpp::NumericVector& vec2,
                            int k = 3, double base = 10,
                            bool NA_rm = false){
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> v2 = Rcpp::as<std::vector<double>>(vec2);
  std::vector<std::vector<double>> cppMat(v1.size(), std::vector<double>(2));
  for (size_t i = 0; i < v1.size(); ++i) {
    cppMat[i][0] = v1[i];
    cppMat[i][1] = v2[i];
  }

  return CppJoinEntropy_Cont(cppMat,static_cast<size_t>(std::abs(k)),base,NA_rm);
};

// [[Rcpp::export]]
double RcppMutualInformation_Cont(const Rcpp::NumericVector& vec1,
                                  const Rcpp::NumericVector& vec2,
                                  int k = 3, int alg = 1,
                                  bool normalize = false,
                                  bool NA_rm = false){
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> v2 = Rcpp::as<std::vector<double>>(vec2);
  std::vector<std::vector<double>> cppMat(v1.size(), std::vector<double>(2));
  for (size_t i = 0; i < v1.size(); ++i) {
    cppMat[i][0] = v1[i];
    cppMat[i][1] = v2[i];
  }

  return CppMutualInformation_Cont(cppMat,static_cast<size_t>(std::abs(k)),alg,normalize,NA_rm);
}

// [[Rcpp::export]]
double RcppConditionalEntropy_Cont(const Rcpp::NumericVector& vec1,
                                   const Rcpp::NumericVector& vec2,
                                   int k = 3, double base = 10,
                                   bool NA_rm = false){
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> v2 = Rcpp::as<std::vector<double>>(vec2);
  return CppConditionalEntropy_Cont(v1,v2,static_cast<size_t>(std::abs(k)),base,NA_rm);
}

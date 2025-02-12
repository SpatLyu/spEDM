#ifndef CppStats_H
#define CppStats_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits
// #include <Rcpp.h>
#include <RcppArmadillo.h>

bool isNA(double value);

bool checkIntNA(int value);

bool checkOneDimVectorHasNaN(const std::vector<double>& vec);

double CppMean(const std::vector<double>& vec,
               bool NA_rm = false);

double CppSum(const std::vector<double>& vec,
              bool NA_rm = false);

double CppMAE(const std::vector<double>& x1,
              const std::vector<double>& x2,
              bool NA_rm = false);

double CppRMSE(const std::vector<double>& x1,
               const std::vector<double>& x2,
               bool NA_rm = false);

std::vector<double> CppAbsDiff(const std::vector<double>& vec1,
                               const std::vector<double>& vec2);

std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false);

double CppVariance(const std::vector<double>& vec, bool NA_rm = false);

double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2,
                     bool NA_rm = false);

double CppDistance(const std::vector<double>& vec1,
                   const std::vector<double>& vec2,
                   bool L1norm = false,
                   bool NA_rm = false);

double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false);

double PartialCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  const std::vector<std::vector<double>>& controls,
                  bool NA_rm = false,
                  bool linear = false);

double PartialCorTrivar(const std::vector<double>& y,
                        const std::vector<double>& y_hat,
                        const std::vector<double>& control,
                        bool NA_rm = false,
                        bool linear = false);

double CppCorSignificance(double r, int n, int k = 0);

std::vector<double> CppCorConfidence(double r, int n, int k = 0,
                                     double level = 0.05);

std::vector<std::size_t> CppKNNIndice(
    const std::vector<std::vector<double>>& embedding_space,
    std::size_t target_idx,
    std::size_t k);

std::vector<std::vector<std::vector<double>>> CppSVD(const std::vector<std::vector<double>>& X);

std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false);

#endif // CppStats_H

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

std::vector<double> CppAbs(const std::vector<double>& vec1,
                           const std::vector<double>& vec2);

std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false);

double CppVariance(const std::vector<double>& vec, bool NA_rm = false);

double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2,
                     bool NA_rm = false);

double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false);

double ArmaPearsonCor(const std::vector<double>& y,
                      const std::vector<double>& y_hat,
                      bool NA_rm = false);

double CppSignificance(double r, int n);

std::vector<double> CppConfidence(double r, int n,
                                  double level = 0.05);

std::vector<std::vector<std::vector<double>>> CppSVD(const std::vector<std::vector<double>>& X);

std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false);

std::vector<double> ArmaLinearTrendRM(const std::vector<double>& vec,
                                      const std::vector<double>& xcoord,
                                      const std::vector<double>& ycoord,
                                      bool NA_rm = false);

#endif // CppStats_H

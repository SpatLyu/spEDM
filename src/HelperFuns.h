#ifndef HelperFuns_H
#define HelperFuns_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits

bool isNA(double value);

bool checkIntNA(int value);

double CppMean(const std::vector<double>& vec,
               bool NA_rm = false);

double CppSum(const std::vector<double>& vec,
              bool NA_rm = false);

std::vector<double> CppAbs(const std::vector<double>& vec1,
                           const std::vector<double>& vec2);

std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false);

#endif // HelperFuns_H

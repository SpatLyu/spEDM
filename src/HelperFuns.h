#ifndef HelperFuns_H
#define HelperFuns_H

#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits

bool isNA(double value);

bool checkIntNA(int value);

double CppMean(const std::vector<double>& vec,
               bool NA_rm = false);

std::vector<double> CppAbs(const std::vector<double>& vec1,
                           const std::vector<double>& vec2);

#endif // HelperFuns_H

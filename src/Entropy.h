#ifndef Entropy_H
#define Entropy_H

#include <cmath>
#include <vector>
#include <numeric>
#include <limits>
#include <map>
#include "CppStats.h"

double CppEntropy_Cont(const std::vector<double>& vec, size_t k,
                       double base = 10, bool NA_rm = false);

double CppJoinEntropy_Cont(const std::vector<std::vector<double>>& mat, size_t k,
                           double base = 10, bool NA_rm = false);

double CppMutualInformation_Cont(const std::vector<std::vector<double>>& mat, size_t k, int alg = 1,
                                 bool normalize = true, bool NA_rm = false);

double CppConditionalEntropy_Cont(const std::vector<double>& vecx,
                                  const std::vector<double>& vecy,
                                  size_t k, double base = 10,
                                  bool NA_rm = false);

double CppEntropy_Disc(const std::vector<double>& vec,
                       double base = 10, bool NA_rm = false);

double CppJoinEntropy_Disc(const std::vector<std::vector<double>>& mat,
                           double base = 10, bool NA_rm = false);

double CppMutualInformation_Disc(const std::vector<std::vector<double>>& mat,
                                 double base = 10, bool NA_rm = false);

double CppConditionalEntropy_Disc(const std::vector<double>& vecx,
                                  const std::vector<double>& vecy,
                                  double base = 10,
                                  bool NA_rm = false);

#endif // Entropy_H

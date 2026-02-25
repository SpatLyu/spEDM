#ifndef Entropy_H
#define Entropy_H

#include <cmath>
#include <vector>
#include <numeric>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "CppStats.h"
#include "CppDistances.h"

double CppEntropy_Cont(const std::vector<double>& vec, size_t k,
                       double base = 10, bool na_rm = true);

double CppJoinEntropy_Cont(const std::vector<std::vector<double>>& mat,
                           const std::vector<int>& columns, size_t k,
                           double base = 10, bool na_rm = true);

double CppMutualInformation_Cont(const std::vector<std::vector<double>>& mat,
                                 const std::vector<int>& columns1,
                                 const std::vector<int>& columns2,
                                 size_t k, int alg = 1,
                                 bool normalize = false, bool na_rm = true);

double CppConditionalEntropy_Cont(const std::vector<std::vector<double>>& mat,
                                  const std::vector<int>& target_columns,
                                  const std::vector<int>& conditional_columns,
                                  size_t k, double base = 10, bool na_rm = true);

double CppEntropy_Disc(const std::vector<double>& vec,
                       double base = 10, bool na_rm = true);

double CppJoinEntropy_Disc(const std::vector<std::vector<double>>& mat,
                           const std::vector<int>& columns,
                           double base = 10, bool na_rm = true);

double CppMutualInformation_Disc(const std::vector<std::vector<double>>& mat,
                                 const std::vector<int>& columns1,
                                 const std::vector<int>& columns2,
                                 double base = 10, bool na_rm = true);

double CppConditionalEntropy_Disc(const std::vector<std::vector<double>>& mat,
                                  const std::vector<int>& target_columns,
                                  const std::vector<int>& conditional_columns,
                                  double base = 10, bool na_rm = true);

#endif // Entropy_H

#ifndef Entropy_H
#define Entropy_H

#include <cmath>
#include <vector>
#include "CppStats.h"

double CppEntropy(const std::vector<double>& vec, size_t k,
                  double base = 10, bool NA_rm = false);

double CppJoinEntropy(const std::vector<std::vector<double>>& mat, size_t k,
                      double base = 10, bool NA_rm = false);

double CppMutualInformation(const std::vector<std::vector<double>>& mat, size_t k, int alg = 1,
                            bool normalize = true, bool NA_rm = false);

double CppConditionalEntropy(const std::vector<double>& vecx,
                             const std::vector<double>& vecy,
                             size_t k, double base = 10,
                             bool NA_rm = false);

#endif // Entropy_H

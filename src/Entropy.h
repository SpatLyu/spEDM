#ifndef Entropy_H
#define Entropy_H

#include <cmath>
#include <vector>
#include "CppStats.h"

double CppEntropy(const std::vector<double>& vec,
                  size_t k, double base = 10,
                  bool L1norm = false, bool NA_rm = false);

double CppJoinEntropy(const std::vector<std::vector<double>>& mat,
                      size_t k, double base = 10,
                      bool L1norm = false, bool NA_rm = false);

#endif // Entropy_H

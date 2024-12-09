#ifndef CppGCCM_H
#define CppGCCM_H

#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>
#include "CppStats.h"
#include "CppUtils.h"
#include <RcppThread.h>

std::vector<std::vector<double>> GCCMLattice(const std::vector<double>& y,
                                             const std::vector<double>& x,
                                             const std::vector<std::vector<int>>& nbmat,
                                             const std::vector<int>& libsizes,
                                             int E);

#endif // CppGCCM_H

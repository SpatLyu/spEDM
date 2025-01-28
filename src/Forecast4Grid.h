#ifndef Forecast4Grid_H
#define Forecast4Grid_H

#include <vector>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& mat,
                                              const std::vector<bool>& lib_indices,
                                              const std::vector<bool>& pred_indices,
                                              const std::vector<int>& E,
                                              double b,
                                              int threads,
                                              bool includeself);

#endif // Forecast4Grid_H

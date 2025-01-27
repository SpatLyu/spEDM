#ifndef Forecast4Lattice_H
#define Forecast4Lattice_H

#include <vector>
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

std::vector<std::vector<double>> Simplex4Lattice(const std::vector<double>& vec,
                                                 const std::vector<std::vector<int>>& nb_vec,
                                                 const std::vector<bool>& lib_indices,
                                                 const std::vector<bool>& pred_indices,
                                                 const std::vector<int>& E,
                                                 double b,
                                                 int threads);

#endif // Forecast4Lattice_H

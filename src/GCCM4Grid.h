#ifndef GCCM4Grid_H
#define GCCM4Grid_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppGridUtils.h"
#include <RcppThread.h>

// Function to perform Simplex Projection Grid
double SimplexProjectionGrid(const std::vector<std::vector<std::vector<double>>>& embeddings,
                             const std::vector<double>& target,
                             const std::vector<bool>& lib_indices,
                             const std::vector<bool>& pred_indices,
                             int num_neighbors);

// GCCMSingle4Grid function
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b);

// GCCM4Grid function
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<int>& lib_sizes,
    const std::vector<std::pair<int, int>>& pred,
    int E,
    int tau = 1,
    int b = 0);

#endif // GCCM4Grid_H

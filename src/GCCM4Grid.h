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
#include "SimplexProjection.h"
#include <RcppThread.h>

// GCCMSingle4Grid function
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b);

// GCCM4Grid function
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix, // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix, // Two dimension matrix of Y variable
    const std::vector<int>& lib_sizes,               // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& pred,    // Indices of spatial units to be predicted
    int E,                                           // Number of dimensions for the attractor reconstruction
    int tau = 1,                                     // Step of spatial lags
    int b = 0,                                       // Number of nearest neighbors to use for prediction
    bool progressbar = true                           // Whether to print the progress bar
);

#endif // GCCM4Grid_H

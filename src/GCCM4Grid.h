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
#include "SMap.h"
#include <RcppThread.h>

/**
 * Calculate the one-dimensional index for a specified position in a 2D grid.
 *
 * This function converts a row and column position in a 2D grid into a corresponding
 * one-dimensional index, assuming the grid is stored in row-major order (i.e., rows are stored sequentially).
 *
 * @param curRow   The current row number (1-based indexing).
 * @param curCol   The current column number (1-based indexing).
 * @param totalRow The total number of rows in the grid.
 * @param totalCol The total number of columns in the grid.
 * @return         The calculated one-dimensional index (0-based indexing).
 */
int LocateGridIndices(int curRow, int curCol,
                      int totalRow, int totalCol);

std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta
);

std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix, // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix, // Two dimension matrix of Y variable
    const std::vector<int>& lib_sizes,               // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& pred,    // Indices of spatial units to be predicted
    int E,                                           // Number of dimensions for the attractor reconstruction
    int tau,                                         // Step of spatial lags
    int b,                                           // Number of nearest neighbors to use for prediction
    bool simplex,                                    // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                    // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                     // Number of threads used from the global pool
    bool progressbar                                 // Whether to print the progress bar
);

#endif // GCCM4Grid_H

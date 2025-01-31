#ifndef SCPCM4Grid_H
#define SCPCM4Grid_H

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
#include "spEDMDataStruct.h"
#include <RcppThread.h>

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size.
 *
 * This function calculates the partial cross mapping between a predictor variable (xEmbedings) and a response
 * variable (yPred) over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings   A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
 * @param yPred        A 1D vector of the response variable's values (spatial cross-section data).
 * @param controls     A 2D matrix that stores the control variables.
 * @param lib_size     The size of the library (number of spatial units) used for prediction.
 * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param conEs        Number of dimensions for the attractor reconstruction with control variables
 * @param totalRow     The total number of rows in the 2D grid.
 * @param totalCol     The total number of columns in the 2D grid.
 * @param b            The number of nearest neighbors to use for prediction.
 * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param cumulate     Whether to cumulate the partial correlations.
 * @param includeself  Whether to include the current state when constructing the embedding vector
 * @return             A vector contains the library size and the corresponding cross mapping and partial cross mapping result.
 */
std::vector<PartialCorRes> SCPCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    int lib_size,
    const std::vector<std::pair<int, int>>& pred,
    const std::vector<int>& conEs,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    bool cumulate,
    bool includeself);

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for multiple library sizes.
 *
 * This function calculates the partial cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
 * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
 *
 * @param xMatrix      A 2D matrix of the predictor variable's values (spatial cross-section data).
 * @param yMatrix      A 2D matrix of the response variable's values (spatial cross-section data).
 * @param zMatrixs     A 2D matrix that stores the control variables.
 * @param lib_sizes    A vector of library sizes (number of spatial units) to use for prediction.
 * @param pred         A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param Es           Number of dimensions for the attractor reconstruction with the x and control variables.
 * @param tau          The step of spatial lags for prediction.
 * @param b            The number of nearest neighbors to use for prediction.
 * @param simplex      If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta        The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads      The number of threads to use for parallel processing.
 * @param cumulate     Whether to cumulate the partial correlations.
 * @param includeself  Whether to include the current state when constructing the embedding vector.
 * @param progressbar  If true, display a progress bar during computation.
 * @return             A 2D vector where each row contains the library size, mean cross mapping result,
 *                     significance, and confidence interval bounds.
 */
std::vector<std::vector<double>> SCPCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,     // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix,     // Two dimension matrix of Y variable
    const std::vector<std::vector<double>>& zMatrixs,    // 2D matrix that stores the control variables
    const std::vector<int>& lib_sizes,                   // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& pred,        // Indices of spatial units to be predicted
    const std::vector<int>& Es,                          // Number of dimensions for the attractor reconstruction with the x and control variables
    int tau,                                             // Step of spatial lags
    int b,                                               // Number of nearest neighbors to use for prediction
    bool simplex,                                        // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                        // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                         // Number of threads used from the global pool
    bool cumulate,                                       // Whether to cumulate the partial correlations
    bool includeself,                                    // Whether to include the current state when constructing the embedding vector
    bool progressbar                                     // Whether to print the progress bar
);

#endif // SCPCM4Grid_H

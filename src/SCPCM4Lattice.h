#ifndef SCPCM4Lattice_H
#define SCPCM4Lattice_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppLatticeUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "spEDMDataStruct.h"
#include <RcppThread.h>

std::vector<PartialCorRes> SCPCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                       // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,        // Neighbor indices vector of the spatial units
    const std::vector<bool>& lib_indices,               // Vector of T/F values (which states to include when searching for neighbors)
    int lib_size,                                       // Size of the library
    int max_lib_size,                                   // Maximum size of the library
    const std::vector<int>& possible_lib_indices,       // Indices of possible library states
    const std::vector<bool>& pred_indices,              // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                      // Number of dimensions for the attractor reconstruction with control variables
    const std::vector<int>& taus,                       // Spatial lag step for constructing lagged state-space vectors with control variables
    int b,                                              // Number of neighbors to use for simplex projection
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    bool cumulate                                       // Whether to cumulate the partial correlations
);

std::vector<std::vector<double>> SCPCM4Lattice(
    const std::vector<double>& x,                       // Spatial cross-section series to cross map from
    const std::vector<double>& y,                       // Spatial cross-section series to cross map to
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**stored by row**)
    const std::vector<std::vector<int>>& nb_vec,        // Neighbor indices vector of the spatial units
    const std::vector<int>& lib_sizes,                  // Vector of library sizes to use
    const std::vector<int>& lib,                        // Vector specifying the library indices
    const std::vector<int>& pred,                       // Vector specifying the prediction indices
    const std::vector<int>& Es,                         // Number of dimensions for the attractor reconstruction with the x and control variables
    const std::vector<int>& taus,                       // Spatial lag step for constructing lagged state-space vectors with the x and control variables
    int b,                                              // Number of nearest neighbors to use for prediction
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                        // Number of threads used from the global pool
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool progressbar                                    // Whether to print the progress bar
);

#endif // SCPCM4Lattice_H

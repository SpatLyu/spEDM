#ifndef GCCM4Lattice_H
#define GCCM4Lattice_H

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
#include <RcppThread.h>

// Function to compute GCCMSingle4Lattice
std::vector<std::pair<int, double>> GCCMSingle4Lattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                      // Time series to be used as the target (should line up with vectors)
    const std::vector<bool>& lib_indices,              // Vector of T/F values (which states to include when searching for neighbors)
    int lib_size,                                      // Size of the library
    int max_lib_size,                                  // Maximum size of the library
    const std::vector<int>& possible_lib_indices,      // Indices of possible library states
    const std::vector<bool>& pred_indices,             // Vector of T/F values (which states to predict from)
    int b                                              // Number of neighbors to use for simplex projection
);

// Function to compute GCCM4Lattice
std::vector<std::vector<double>> GCCM4Lattice(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                       // Spatial cross-section series to cross map to
    const std::vector<int>& lib_sizes,                  // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& lib,        // Matrix (n x 2) using n sequences of data to construct libraries
    const std::vector<std::pair<int, int>>& pred,       // Matrix (n x 2) using n sequences of data to predict from
    int E,                                              // Number of dimensions for the attractor reconstruction
    int tau = 1,                                        // Spatial lag for the lagged-vector construction
    int b = 0,                                          // Number of nearest neighbors to use for prediction
    bool progressbar = true                             // Whether to print the progress bar
);

#endif // GCCM4Lattice_H

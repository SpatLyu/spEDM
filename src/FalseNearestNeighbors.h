#ifndef FalseNearestNeighbors_H
#define FalseNearestNeighbors_H

#include <vector>
#include <cmath>
#include <limits>
#include "CppStats.h"
#include <RcppThread.h>

/*
 * Compute the False Nearest Neighbors (FNN) ratio for spatial cross-sectional data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - Rtol: Relative threshold (default 10.0). If the change in the added dimension is
 *         large relative to the E1 distance, the neighbor is considered false.
 * - Atol: Absolute threshold (default 2.0). If the added dimension changes too much
 *         absolutely, the neighbor is also considered false.
 * - L1norm: Whether to use Manhattan (L1) distance instead of Euclidean (L2).
 *
 * Returns:
 * - A double value indicating the proportion of false nearest neighbors (0–1).
 *   If no valid pairs are found, returns NaN.
 */
double CppSingleFNN(const std::vector<std::vector<double>>& embedding,
                    size_t E1,
                    size_t E2,
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false);

#endif // FalseNearestNeighbors_H

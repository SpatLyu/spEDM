#ifndef SignatureProjection_H
#define SignatureProjection_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "NumericUtils.h"
#include <RcppThread.h>
/**
 * @brief Predicts signature vectors for a subset of target points using weighted nearest neighbors.
 *
 * This function performs local weighted prediction in the signature space as follows:
 *   1. For each prediction index `p` in `pred_indices`, find its `num_neighbors` nearest neighbors
 *      among `lib_indices` based on distances in `Dx[p][*]`, ignoring NaN distances.
 *   2. Compute exponential weights scaled by the total distance sum to emphasize close points.
 *      If all distances are zero, uniform weights are used instead.
 *   3. For each dimension of the signature space:
 *        - Count how many neighbor signatures are exactly zero.
 *        - If the zero count exceeds `zero_tolerance`, set the predicted value to 0.
 *        - Otherwise, compute a weighted average of valid (non-NaN) neighbor signatures.
 *   4. Predictions are stored and updated only for indices in `pred_indices`; other entries remain undefined (NaN).
 *
 * Parallelization:
 *   - Controlled by the parameter `threads`.
 *   - If `threads <= 1`, computation is serial (standard for-loop).
 *   - Otherwise, the loop over prediction indices is executed in parallel via RcppThread::parallelFor.
 *
 * @param SMy            Signature space of the target variable Y. Shape: (N_obs, E−1)
 * @param Dx             Distance matrix from prediction points to library points. Shape: (SMy.size(), SMy.size())
 * @param lib_indices    Indices of valid library points used for neighbor search (subset of [0, SMy.size())).
 * @param pred_indices   Indices of points to predict (subset of [0, SMy.size())).
 * @param num_neighbors  Number of nearest neighbors to use. If <= 0, defaults to E+1.
 * @param zero_tolerance Maximum allowed zero values per dimension before forcing prediction to zero.
 *                       If <= 0, defaults to E−1.
 * @param threads        Number of threads to use. If <= 1, runs serially; otherwise runs parallel.
 *
 * @return A matrix of predicted signature vectors, sized SMy.size() × (E−1).
 */
std::vector<std::vector<double>> SignatureProjection(
    const std::vector<std::vector<double>>& SMy,
    const std::vector<std::vector<double>>& Dx,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    int num_neighbors = 0,   /* = std::numeric_limits<int>::min() */
    int zero_tolerance = 0,  /* = std::numeric_limits<int>::max() */
    size_t threads = 1
);

#endif // SignatureProjection_H

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

// [[Rcpp::depends(RcppThread)]]

/**
 * @brief Computes the partial correlation between a spatial cross-section series and its prediction
 *        using the Simplex Projection method, incorporating control variables in a grid-based spatial setting.
 *
 * This function reconstructs the state-space and applies Simplex Projection prediction while accounting
 * for control variables in a spatial grid data. The process can be cumulative or independent in incorporating
 * control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Spatial cross-section series used as the prediction target.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param lib_indices: Boolean vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Boolean vector indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for Simplex Projection.
 * @param nrow: Number of rows in the input spatial grid data.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSimplex4Grid(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2,std::numeric_limits<double>::quiet_NaN());

  if (cumulate) {
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      } else {
        temp_pred = SimplexProjectionPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i]);
      }
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls]);
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
    }
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],taus[i]);
      temp_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1]);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCor(target,target_pred,con_pred,true,false);
    }
  }

  return rho;
}

/**
 * @brief Computes the partial correlation between a spatial cross-section series and its prediction
 *        using the S-Map method, incorporating control variables in a grid-based spatial setting.
 *
 * This function reconstructs the state-space and applies S-Map prediction while accounting for
 * control variables in a spatial grid data. The process can be cumulative or independent in incorporating
 * control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Spatial cross-section series used as the prediction target.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param lib_indices: Boolean vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Boolean vector indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for Simplex Projection.
 * @param nrow: Number of rows in the input spatial grid data.
 * @param theta: Weighting parameter for distances in the S-Map method.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSMap4Grid(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<bool>& lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    double theta,
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2,std::numeric_limits<double>::quiet_NaN());

  if (cumulate){
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      } else {
        temp_pred = SMapPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i], theta);
      }
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls], theta);
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
    }
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],taus[i]);
      temp_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1], theta);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCor(target,target_pred,con_pred,true,false);
    }
  }

  return rho;
}

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size.
 *
 * This function calculates the partial cross mapping between a predictor variable (xEmbedings) and a response
 * variable (yPred) over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings           A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
 * @param yPred                A 1D vector of the response variable's values (spatial cross-section data).
 * @param controls             A 2D matrix that stores the control variables.
 * @param lib_sizes            A vector of two integers, where the first element is the row-wise library size and the second element is the column-wise library size.
 * @param possible_lib_indices A boolean vector indicating which spatial units are valid for inclusion in the library.
 * @param pred_indices         A boolean vector indicating which spatial units to be predicted.
 * @param conEs                Number of dimensions for the attractor reconstruction with control variables
 * @param taus:                Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param totalRow             The total number of rows in the 2D grid.
 * @param totalCol             The total number of columns in the 2D grid.
 * @param b                    The numbers of nearest neighbors to use for prediction.
 * @param simplex              If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta                The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads              The number of threads to use for parallel processing.
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * @param cumulate             Whether to cumulate the partial correlations.
 * @param row_size_mark        If ture, use the row-wise libsize to mark the libsize; if false, use col-wise libsize.
 *
 * @return  A vector contains the library size and the corresponding cross mapping and partial cross mapping result.
 */
std::vector<PartialCorRes> SCPCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& possible_lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate,
    bool row_size_mark)
{
  // Extract row-wise and column-wise library sizes
  const int lib_size_row = lib_sizes[0];
  const int lib_size_col = lib_sizes[1];

  // Determine the marked libsize
  const int libsize = row_size_mark ? lib_size_row : lib_size_col;
  const int half_lib_size = (lib_size_row * lib_size_col) / 2;

  // Precompute valid (r, c) pairs
  std::vector<std::pair<int, int>> valid_indices;
  // valid_indices.reserve((totalRow - lib_size_row + 1) * (totalCol - lib_size_col + 1));
  for (int r = 1; r <= totalRow - lib_size_row + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size_col + 1; ++c) {
      valid_indices.emplace_back(r, c);
    }
  }

  // // Initialize the result container with the same size as valid_indices
  // std::vector<PartialCorRes> x_xmap_y;
  // x_xmap_y.resize(valid_indices.size());

  // Initialize the result container with the same size as valid_indices,
  // and optionally set default values using the constructor of PartialCorRes
  std::vector<PartialCorRes> x_xmap_y(valid_indices.size());

  // Unified processing logic
  auto process = [&](size_t idx) {
    const int r = valid_indices[idx].first;
    const int c = valid_indices[idx].second;

    // Initialize library indices and count the number of the nan value together
    std::vector<bool> lib_indices(totalRow * totalCol, false);
    int na_count = 0;

    for (int ii = r; ii < r + lib_size_row; ++ii) {
      for (int jj = c; jj < c + lib_size_col; ++jj) {
        const int index = (ii - 1) * totalCol + (jj - 1);
        if (possible_lib_indices[index]) {
          lib_indices[index] = true;
          if (std::isnan(yPred[index])) {
            ++na_count;
          }
        }
      }
    }

    // Directly initialize std::vector<double> with two NaN values
    std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());

    if (na_count <= half_lib_size) {
      // Run partial cross map and store results
      if (simplex) {
        rho = PartialSimplex4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, cumulate);
      } else {
        rho = PartialSMap4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, theta, cumulate);
      }
    }

    // Directly assign a PartialCorRes struct with the three values
    x_xmap_y[idx] = PartialCorRes(libsize, rho[0], rho[1]);
  };

  // Parallel coordination
  if (parallel_level == 0) {
    RcppThread::parallelFor(0, valid_indices.size(), process, threads);
  } else {
    for (size_t i = 0; i < valid_indices.size(); ++i) {
      process(i);
    }
  }

  return x_xmap_y;
}

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size.
 *
 * This function follows the same library construction logic as SCPCMSingle4Lattice, where libraries
 * are created by selecting consecutive indices from possible_lib_indices with possible wraparound.
 *
 * @param xEmbedings           State-space embeddings for the predictor variable (each row is a spatial vector)
 * @param yPred                Target spatial cross-section series
 * @param controls             Control variables stored by row
 * @param lib_size             Number of consecutive spatial units to include in each library
 * @param max_lib_size         Maximum possible library size (total valid spatial units)
 * @param possible_lib_indices Integer vector indicating the indices of eligible spatial units for library construction
 * @param pred_indices         Boolean vector indicating spatial units to predict
 * @param conEs                Embedding dimensions for control variables
 * @param taus                 Spatial lag steps for control variable embeddings
 * @param b                    Number of nearest neighbors for prediction
 * @param totalRow             Total rows in spatial grid
 * @param totalCol             Total columns in spatial grid
 * @param simplex              Use simplex projection if true, S-mapping if false
 * @param theta                Distance weighting parameter for S-mapping
 * @param threads              The number of threads to use for parallel processing
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`
 * @param cumulate             Enable cumulative partial correlations
 *
 * @return Vector of PartialCorRes containing mapping results for each library configuration
 */
std::vector<PartialCorRes> SCPCMSingle4GridOneDim(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    int lib_size,
    int max_lib_size,
    const std::vector<int>& possible_lib_indices,
    const std::vector<bool>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate) {
  int n = yPred.size();

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<bool> lib_indices(n, false);
    for (int idx : possible_lib_indices) {
      lib_indices[idx] = true;
    }

    // Check if more than half of the library is NA
    int na_count = 0;
    for (size_t i = 0; i < lib_indices.size(); ++i) {
      if (lib_indices[i] && std::isnan(yPred[i])) {
        ++na_count;
      }
    }

    std::vector<PartialCorRes> x_xmap_y;
    std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());
    if (na_count <= max_lib_size / 2.0) {
      // Run partial cross map and store results
      if (simplex) {
        rho = PartialSimplex4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, cumulate);
      } else {
        rho = PartialSMap4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, theta, cumulate);
      }
    }

    x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    return x_xmap_y;
  } else if (parallel_level == 0) {
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(i);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(i);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(i);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<PartialCorRes> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      std::vector<bool> lib_indices(n, false);
      std::vector<int> local_lib_indices = valid_lib_indices[i];
      for(int& li : local_lib_indices){
        lib_indices[possible_lib_indices[li]] = true;
      }

      // Check if more than half of the library is NA
      int na_count = 0;
      for (size_t ii = 0; ii < lib_indices.size(); ++ii) {
        if (lib_indices[ii] && std::isnan(yPred[ii])) {
          ++na_count;
        }
      }

      std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());
      if (na_count <= max_lib_size / 2.0) {
        // Run partial cross map and store results
        if (simplex) {
          rho = PartialSimplex4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, cumulate);
        } else {
          rho = PartialSMap4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, theta, cumulate);
        }
      }

      // Directly initialize a PartialCorRes struct with the three values
      PartialCorRes result(lib_size, rho[0], rho[1]);
      x_xmap_y[i] = result;
    }, threads);

    return x_xmap_y;
  } else {
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(i);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(i);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(i);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<PartialCorRes> x_xmap_y(valid_lib_indices.size());

    // Iterate through precomputed valid_lib_indices
    for (size_t i = 0; i < valid_lib_indices.size(); ++i){
      std::vector<bool> lib_indices(n, false);
      std::vector<int> local_lib_indices = valid_lib_indices[i];
      for(int& li : local_lib_indices){
        lib_indices[possible_lib_indices[li]] = true;
      }

      // Check if more than half of the library is NA
      int na_count = 0;
      for (size_t ii = 0; ii < lib_indices.size(); ++ii) {
        if (lib_indices[ii] && std::isnan(yPred[ii])) {
          ++na_count;
        }
      }

      std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());
      if (na_count <= max_lib_size / 2.0) {
        // Run partial cross map and store results
        if (simplex) {
          rho = PartialSimplex4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, cumulate);
        } else {
          rho = PartialSMap4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, taus, b, totalRow, theta, cumulate);
        }
      }

      // Directly initialize a PartialCorRes struct with the three values
      PartialCorRes result(lib_size, rho[0], rho[1]);
      x_xmap_y[i] = result;
    }

    return x_xmap_y;
  }
}

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for multiple library sizes.
 *
 * This function estimates the cross mapping and partial cross mapping between predictor variables (`xMatrix`) and response
 * variables (`yMatrix`) over a 2D spatial grid, incorporating control variables (`zMatrixs`). It supports both Simplex Projection
 * and S-Mapping, with options for parallel computation and progress tracking.
 *
 * Parameters:
 * - xMatrix: A 2D matrix of predictor variable values (spatial cross-section data).
 * - yMatrix: A 2D matrix of response variable values (spatial cross-section data).
 * - zMatrixs: A 2D matrix storing the control variables.
 * - lib_sizes: A 2D vector where the first sub-vector contains row-wise library sizes and the second sub-vector contains column-wise library sizes.
 * - lib: A vector of pairs representing the indices (row, column) of spatial units to be the library.
 * - pred: A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * - Es: A vector specifying the embedding dimensions for attractor reconstruction using `xMatrix` and control variables.
 * - taus: A vector specifying the spatial lag steps for constructing lagged state-space vectors with control variables.
 * - b: A vector specifying the numbers of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use Simplex Projection (true) or S-Mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-Mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * - cumulate: Boolean flag indicating whether to cumulate partial correlations.
 * - progressbar: Boolean flag indicating whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean pearson cross-mapping correlation.
 *      - The statistical significance of the pearson correlation.
 *      - The upper bound of the pearson correlation confidence interval.
 *      - The lower bound of the pearson correlation confidence interval.
 *      - The mean partial cross-mapping partial correlation.
 *      - The statistical significance of the partial correlation.
 *      - The upper bound of the partial correlation confidence interval.
 *      - The lower bound of the partial correlation confidence interval.
 */
std::vector<std::vector<double>> SCPCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,     // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix,     // Two dimension matrix of Y variable
    const std::vector<std::vector<double>>& zMatrixs,    // 2D matrix that stores the control variables
    const std::vector<std::vector<int>>& lib_sizes,      // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& lib,         // Indices of spatial units to be the library
    const std::vector<std::pair<int, int>>& pred,        // Indices of spatial units to be predicted
    const std::vector<int>& Es,                          // Number of dimensions for the attractor reconstruction with the x and control variables
    const std::vector<int>& taus,                        // Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
    const std::vector<int>& b,                           // Numbers of nearest neighbors to use for prediction
    bool simplex,                                        // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                        // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                         // Number of threads used from the global pool
    int parallel_level,                                  // Level of parallel computing: 0 for `lower`, 1 for `higher`
    bool cumulate,                                       // Whether to cumulate the partial correlations
    bool progressbar                                     // Whether to print the progress bar
) {
  // If b is not provided correctly, default it to E + 2
  std::vector<int> bs = b;
  for (size_t i = 0; i < bs.size(); ++i){
    if (bs[i] <= 0) {
      bs[i] = Es[i] + 2;
    }
  }

  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  int taux = taus[0];
  std::vector<int> contaus = taus;
  contaus.erase(contaus.begin());

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, Ex, taux);

  int n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = zMatrixs.size();
  }

  // Ensure the maximum value does not exceed totalRow or totalCol
  int max_lib_size_row = totalRow;
  int max_lib_size_col = totalCol;

  // Extract row-wise and column-wise library sizes
  std::vector<int> row_lib_sizes = lib_sizes[0];
  std::vector<int> col_lib_sizes = lib_sizes[1];

  // Transform to ensure no size exceeds max_lib_size_row or max_lib_size_col
  std::transform(row_lib_sizes.begin(), row_lib_sizes.end(), row_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_row); });
  std::transform(col_lib_sizes.begin(), col_lib_sizes.end(), col_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_col); });

  // Remove duplicates in row-wise and column-wise library sizes
  std::sort(row_lib_sizes.begin(), row_lib_sizes.end());
  row_lib_sizes.erase(std::unique(row_lib_sizes.begin(), row_lib_sizes.end()), row_lib_sizes.end());
  std::sort(col_lib_sizes.begin(), col_lib_sizes.end());
  col_lib_sizes.erase(std::unique(col_lib_sizes.begin(), col_lib_sizes.end()), col_lib_sizes.end());

  // // Generate unique pairs of row-wise and column-wise library sizes
  // std::vector<std::pair<int, int>> unique_lib_size_pairs;
  // for (int row_size : row_lib_sizes) {
  //   for (int col_size : col_lib_sizes) {
  //     unique_lib_size_pairs.emplace_back(row_size, col_size);
  //   }
  // }

  // Generate unique pairs of row-wise and column-wise library sizes
  std::vector<std::pair<int, int>> unique_lib_size_pairs;

  // Determine which library size vector is longer
  int row_size_count = row_lib_sizes.size();
  int col_size_count = col_lib_sizes.size();
  int min_size = std::min(row_size_count, col_size_count);
  // int max_size = std::max(row_size_count, col_size_count);

  // Fill unique_lib_size_pairs based on the shorter vector
  for (int i = 0; i < min_size; ++i) {
    unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes[i]);
  }

  bool row_size_mark = true;
  // Handle the excess elements for the longer vector
  if (row_size_count > col_size_count) {
    for (int i = min_size; i < row_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes.back()); // Pair with the max value of col_lib_sizes
    }
  }

  if (row_size_count < col_size_count) {
    for (int i = min_size; i < col_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes.back(), col_lib_sizes[i]); // Pair with the max value of row_lib_sizes
    }
    row_size_mark = false;
  }

  // Set library indices
  std::vector<bool> lib_indices(totalRow * totalCol, false);
  for (const auto& l : lib) {
    lib_indices[LocateGridIndices(l.first + 1, l.second + 1, totalRow, totalCol)] = true;
  }

  // Set prediction indices
  std::vector<bool> pred_indices(totalRow * totalCol, false);
  for (const auto& p : pred) {
    pred_indices[LocateGridIndices(p.first + 1, p.second + 1, totalRow, totalCol)] = true;
  }

  // Exclude NA values in yPred from the library and prediction indices
  for (size_t i = 0; i < yPred.size(); ++i) {
    if (std::isnan(yPred[i])) {
      lib_indices[i] = false;
      pred_indices[i] = false;
    }
  }

  // Local results for each library
  std::vector<std::vector<PartialCorRes>> local_results(unique_lib_size_pairs.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
      for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        local_results[i] = SCPCMSingle4Grid(
          xEmbedings,
          yPred,
          zMatrixs,
          {lib_size_row, lib_size_col},
          lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate,
          row_size_mark);
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        local_results[i] = SCPCMSingle4Grid(
          xEmbedings,
          yPred,
          zMatrixs,
          {lib_size_row, lib_size_col},
          lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate,
          row_size_mark);
      }
    }
  } else {
    // Perform the operations using RcppThread
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
      RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        local_results[i] = SCPCMSingle4Grid(
          xEmbedings,
          yPred,
          zMatrixs,
          {lib_size_row, lib_size_col},
          lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate,
          row_size_mark);
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
      RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        local_results[i] = SCPCMSingle4Grid(
          xEmbedings,
          yPred,
          zMatrixs,
          {lib_size_row, lib_size_col},
          lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate,
          row_size_mark);
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<PartialCorRes> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int and store second and third values as pairs
  std::map<int, std::vector<std::pair<double, double>>> grouped_results;

  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].emplace_back(result.second, result.third);
  }

  std::vector<std::vector<double>> final_results;

  // Compute the mean of second and third values for each group
  for (const auto& group : grouped_results) {
    std::vector<double> second_values, third_values;

    for (const auto& val : group.second) {
      second_values.push_back(val.first);
      third_values.push_back(val.second);
    }

    double mean_second = CppMean(second_values, true);
    double mean_third = CppMean(third_values, true);

    final_results.push_back({static_cast<double>(group.first), mean_second, mean_third});
  }

  int n = pred.size();
  // Compute significance and confidence intervals for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho_second = final_results[i][1];
    double rho_third = final_results[i][2];

    // Compute significance and confidence interval for second value
    double significance_second = CppCorSignificance(rho_second, n);
    std::vector<double> confidence_interval_second = CppCorConfidence(rho_second, n);

    // Compute significance and confidence interval for third value
    double significance_third = CppCorSignificance(rho_third, n, n_confounds);
    std::vector<double> confidence_interval_third = CppCorConfidence(rho_third, n, n_confounds);

    // Append computed statistical values to the result
    final_results[i].push_back(significance_second);
    final_results[i].push_back(confidence_interval_second[0]);
    final_results[i].push_back(confidence_interval_second[1]);

    final_results[i].push_back(significance_third);
    final_results[i].push_back(confidence_interval_third[0]);
    final_results[i].push_back(confidence_interval_third[1]);
  }

  return final_results;
}

/**
 * Performs Grid-based Spatially Convergent Partial Cross Mapping (SCPCM).
 *
 * Parameters:
 * - xMatrix: 2D grid of predictor variable values (row-major order)
 * - yMatrix: 2D grid of target variable values (row-major order)
 * - zMatrixs: Control variables stored as 2D grids (vector of vectors)
 * - lib_size: Number of consecutive spatial units to include in each library
 * - lib: Vector specifying library location as (row, col) pairs
 * - pred: Vector specifying prediction locations
 * - Es: Embedding dimensions for x and control variables
 * - taus: Spatial lag steps for x and control variables
 * - b: Numbers of nearest neighbors for prediction
 * - simplex: Use simplex projection (true) or S-mapping (false)
 * - theta: Distance weighting parameter for S-mapping
 * - threads: Number of parallel computation threads
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * - cumulate: Enable cumulative partial correlations
 * - progressbar: Display progress bar during computation
 *
 * Returns:
 *   2D vector containing:
 *     - Library size
 *     - Mean cross-map correlation (rho)
 *     - Rho significance
 *     - Rho upper CI
 *     - Rho lower CI
 *     - Mean partial correlation
 *     - Partial correlation significance
 *     - Partial upper CI
 *     - Partial lower CI
 */
std::vector<std::vector<double>> SCPCM4GridOneDim(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<std::vector<double>>& zMatrixs,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    const std::vector<int>& Es,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    bool cumulate,
    bool progressbar
){
  // If b is not provided correctly, default it to E + 2
  std::vector<int> bs = b;
  for (size_t i = 0; i < bs.size(); ++i){
    if (bs[i] <= 0) {
      bs[i] = Es[i] + 2;
    }
  }

  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  int taux = taus[0];
  std::vector<int> contaus = taus;
  contaus.erase(contaus.begin());

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, Ex, taux);

  int n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = zMatrixs.size();
  }

  std::vector<int> possible_lib_indices;
  for (size_t i = 0; i < lib.size(); ++i) {
    int LibIndice = lib[i] - 1;
    if (!std::isnan(yPred[LibIndice])) {
      possible_lib_indices.push_back(LibIndice);
    }
  }
  int max_lib_size = static_cast<int>(possible_lib_indices.size()); // Maximum lib size

  // Initialize pred_indices with all false
  std::vector<bool> pred_indices(totalRow*totalCol, false);
  // Convert pred (1-based in R) to 0-based indices, exclude yPred NA and set corresponding positions to true
  for (size_t i = 0; i < pred.size(); ++i) {
    int PreIndice = pred[i] - 1;
    if (!std::isnan(yPred[PreIndice])) {
      pred_indices[PreIndice] = true;
    }
  }

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is Ex + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, Ex + 2); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Local results for each library
  std::vector<std::vector<PartialCorRes>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = SCPCMSingle4GridOneDim(
          xEmbedings,
          yPred,
          zMatrixs,
          unique_lib_sizes[i],
          max_lib_size,
          possible_lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = SCPCMSingle4GridOneDim(
          xEmbedings,
          yPred,
          zMatrixs,
          unique_lib_sizes[i],
          max_lib_size,
          possible_lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
      }
    }
  } else {
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = SCPCMSingle4GridOneDim(
          xEmbedings,
          yPred,
          zMatrixs,
          lib_size,
          max_lib_size,
          possible_lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = SCPCMSingle4GridOneDim(
          xEmbedings,
          yPred,
          zMatrixs,
          lib_size,
          max_lib_size,
          possible_lib_indices,
          pred_indices,
          conEs,
          contaus,
          bs,
          totalRow,
          totalCol,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<PartialCorRes> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int and store second and third values as pairs
  std::map<int, std::vector<std::pair<double, double>>> grouped_results;

  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].emplace_back(result.second, result.third);
  }

  std::vector<std::vector<double>> final_results;

  // Compute the mean of second and third values for each group
  for (const auto& group : grouped_results) {
    std::vector<double> second_values, third_values;

    for (const auto& val : group.second) {
      second_values.push_back(val.first);
      third_values.push_back(val.second);
    }

    double mean_second = CppMean(second_values, true);
    double mean_third = CppMean(third_values, true);

    final_results.push_back({static_cast<double>(group.first), mean_second, mean_third});
  }

  int n = pred.size();
  // Compute significance and confidence intervals for each result
  for (size_t i = 0; i < final_results.size(); ++i) {
    double rho_second = final_results[i][1];
    double rho_third = final_results[i][2];

    // Compute significance and confidence interval for second value
    double significance_second = CppCorSignificance(rho_second, n);
    std::vector<double> confidence_interval_second = CppCorConfidence(rho_second, n);

    // Compute significance and confidence interval for third value
    double significance_third = CppCorSignificance(rho_third, n, n_confounds);
    std::vector<double> confidence_interval_third = CppCorConfidence(rho_third, n, n_confounds);

    // Append computed statistical values to the result
    final_results[i].push_back(significance_second);
    final_results[i].push_back(confidence_interval_second[0]);
    final_results[i].push_back(confidence_interval_second[1]);

    final_results[i].push_back(significance_third);
    final_results[i].push_back(confidence_interval_third[0]);
    final_results[i].push_back(confidence_interval_third[1]);
  }

  return final_results;
}

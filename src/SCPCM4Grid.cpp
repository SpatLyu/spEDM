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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

std::vector<double> PartialSimplex4Grid(
    const std::vector<std::vector<double>>& vectors,    // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                  // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**each matirx stored by row**)
    const std::vector<bool>& lib_indices,               // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,              // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                      // Number of dimensions for the attractor reconstruction with control variables
    int num_neighbors,                                  // Number of neighbors to use for simplex projection
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool includeself                                    // Whether to include the current state when constructing the embedding vector
){
  int nrow = controls[0].size();
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate) {
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors);
      } else {
        temp_pred = SimplexProjectionPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors);
      }
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],includeself);
    }

    std::vector<double> con_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors);
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors);
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],includeself);
      temp_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCor(target,target_pred,con_pred,true,false);
  }
  return rho;
}

std::vector<double> PartialSMap4Grid(
    const std::vector<std::vector<double>>& vectors,    // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& target,                  // Spatial cross-section series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls,   // Cross-sectional data of control variables (**each variable stored by one row**)
    const std::vector<bool>& lib_indices,               // Vector of T/F values (which states to include when searching for neighbors)
    const std::vector<bool>& pred_indices,              // Vector of T/F values (which states to predict from)
    const std::vector<int>& conEs,                      // Number of dimensions for the attractor reconstruction with control variables
    int num_neighbors,                                  // Number of neighbors to use for simplex projection
    double theta,                                       // Weighting parameter for distances
    bool cumulate,                                      // Whether to cumulate the partial correlations
    bool includeself                                    // Whether to include the current state when constructing the embedding vector
){
  int nrow = controls[0].size();
  int n_controls = controls.size();
  std::vector<double> rho(2);

  if (cumulate){
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      } else {
        temp_pred = SMapPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      }
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],includeself);
    }

    std::vector<double> con_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors, theta);
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
  } else {
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_conmat;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors, theta);
      temp_conmat = GridVec2Mat(temp_pred,nrow);
      temp_embedding = GenGridEmbeddings(temp_conmat,conEs[i],includeself);
      temp_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors, theta);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

    rho[0] = PearsonCor(target,target_pred,true);
    rho[1] = PartialCor(target,target_pred,con_pred,true,false);
  }

  return rho;
}

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
    bool includeself)
{
  std::vector<PartialCorRes> x_xmap_y;
  std::vector<double> rho;

  for (int r = 1; r <= totalRow - lib_size + 1; ++r) {
    for (int c = 1; c <= totalCol - lib_size + 1; ++c) {

      // Initialize prediction and library indices
      std::vector<bool> pred_indices(totalRow * totalCol, false);
      std::vector<bool> lib_indices(totalRow * totalCol, false);

      // Set prediction indices
      for (const auto& p : pred) {
        pred_indices[LocateGridIndices(p.first, p.second, totalRow, totalCol)] = true;
      }

      // Exclude NA values in yPred from prediction indices
      for (size_t i = 0; i < yPred.size(); ++i) {
        if (std::isnan(yPred[i])) {
          pred_indices[i] = false;
        }
      }

      // Set library indices
      for (int i = r; i < r + lib_size; ++i) {
        for (int j = c; j < c + lib_size; ++j) {
          lib_indices[LocateGridIndices(i, j, totalRow, totalCol)] = true;
        }
      }

      // Check if more than half of the library is NA
      int na_count = 0;
      for (size_t i = 0; i < lib_indices.size(); ++i) {
        if (lib_indices[i] && std::isnan(yPred[i])) {
          ++na_count;
        }
      }
      if (na_count > (lib_size * lib_size) / 2) {
        continue;
      }

      // Run partial cross map and store results
      if (simplex) {
        rho = PartialSimplex4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, b, cumulate, includeself);
      } else {
        rho = PartialSMap4Grid(xEmbedings, yPred, controls, lib_indices, pred_indices, conEs, b, theta, cumulate, includeself);
      }
      x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    }
  }

  return x_xmap_y;
}

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
) {
  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  // If b is not provided correctly, default it to Ex + 2
  if (b <= 0) {
    b = Ex + 2;
  }

  size_t threads_sizet = static_cast<size_t>(threads);
  unsigned int max_threads = std::thread::hardware_concurrency();
  threads_sizet = std::min(static_cast<size_t>(max_threads), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> xEmbedings = GenGridEmbeddings(xMatrix, Ex, includeself);

  int n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = zMatrixs.size();
  }

  // Ensure the maximum value does not exceed totalRow or totalCol
  int max_lib_size = std::max(totalRow, totalCol);

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is E + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, E + 2); });

  // Remove duplicates
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Initialize the result container
  std::vector<PartialCorRes> x_xmap_y;

  // // Iterate over each library size
  // for (int lib_size : unique_lib_sizes) {
  //   // Perform single grid partial cross-mapping for the current library size
  //   auto results = SCPCMSingle4Grid(
  //     xEmbedings,
  //     yPred,
  //     zMatrixs,
  //     lib_size,
  //     pred,
  //     conEs,
  //     totalRow,
  //     totalCol,
  //     b,
  //     simplex,
  //     theta,
  //     cumulate,
  //     includeself);
  //
  //   // Append the results to the main result container
  //   x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
  // }

  // Perform the operations using RcppThread
  if (progressbar) {
    RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = SCPCMSingle4Grid(
        xEmbedings,
        yPred,
        zMatrixs,
        lib_size,
        pred,
        conEs,
        totalRow,
        totalCol,
        b,
        simplex,
        theta,
        cumulate,
        includeself);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
      int lib_size = unique_lib_sizes[i];
      auto results = SCPCMSingle4Grid(
        xEmbedings,
        yPred,
        zMatrixs,
        lib_size,
        pred,
        conEs,
        totalRow,
        totalCol,
        b,
        simplex,
        theta,
        cumulate,
        includeself);
      x_xmap_y.insert(x_xmap_y.end(), results.begin(), results.end());
    }, threads_sizet);
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

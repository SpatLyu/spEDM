#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>
#include "CppStats.h"
#include "CppUtils.h"
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

// Function to perform GCCM Lattice
std::vector<std::vector<double>> GCCMLattice(const std::vector<double>& y,
                                             const std::vector<double>& x,
                                             const std::vector<std::vector<int>>& nbmat,
                                             const std::vector<int>& libsizes,
                                             int E) {
  // Check if the maximum and minimum values of libsizes are within the valid range
  int y_size = y.size();

  // Construct a new libsizevec with valid libsize values
  std::vector<int> libsizevec;
  // libsizevec.push_back(y_size - 2);
  for (int libsize : libsizes) {
    if (libsize >= 1 && libsize <= y_size - 2) {
      libsizevec.push_back(libsize);
    }
  }
  // libsizevec.push_back(1);

  // Sort libsizevec in descending order
  std::sort(libsizevec.begin(), libsizevec.end(), std::greater<int>());

  // Initialize the result vector
  std::vector<std::vector<double>> result;

  // Perform the operations for each libsize in libsizevec using RcppThread
  RcppThread::parallelFor(0, libsizevec.size(), [&](size_t i) {
    int libsize = libsizevec[i];

    // Perform Simplex Projection
    std::vector<double> y_hat = SimplexProjection(y, x, nbmat, libsize, E);

    // Calculate the mean of y_hat
    double y_hat_mean = CppMean(y_hat, true);

    // Calculate the Pearson correlation coefficient
    double rho = PearsonCor(y, y_hat, true);

    // Store the results in the result vector
    result.push_back({static_cast<double>(libsize), y_hat_mean, rho});
  });

  // Calculate significance and confidence interval for each result
  for (size_t i = 0; i < result.size(); ++i) {
    double rho = result[i][2];
    double significance = CppSignificance(rho, y_size);
    std::vector<double> confidence_interval = CppConfidence(rho, y_size);

    result[i].push_back(significance);
    result[i].push_back(confidence_interval[0]);
    result[i].push_back(confidence_interval[1]);
  }

  return result;
}

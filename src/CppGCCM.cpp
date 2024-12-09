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

// Function to calculate multiple rho values for each libsize
std::vector<double> CalculateRhoForLibsize(const std::vector<double>& y,
                                           const std::vector<double>& x,
                                           const std::vector<std::vector<int>>& nbmat,
                                           int libsize,
                                           int E) {
  std::vector<double> rhos;
  size_t n = y.size();
  size_t libsize_t = static_cast<size_t>(libsize); // Convert libsize to size_t

  for (size_t start = 0; start < n; ++start) {
    std::vector<double> y_subset;
    std::vector<double> x_subset;
    std::vector<bool> visited(n, false);
    visited[start] = true;

    // Collect neighbors based on nbmat
    std::vector<size_t> current_neighbors = {start};
    while (y_subset.size() < libsize_t && !current_neighbors.empty()) {
      std::vector<size_t> next_neighbors;
      for (size_t neighbor : current_neighbors) {
        for (size_t i = 0; i < n; ++i) {
          if (nbmat[neighbor][i] == 1 && !visited[i]) {
            y_subset.push_back(y[i]);
            x_subset.push_back(x[i]);
            visited[i] = true;
            if (y_subset.size() < libsize_t) {
              next_neighbors.push_back(i);
            }
          }
        }
      }
      current_neighbors = next_neighbors;
    }

    // Ensure the subset size is at least libsize
    if (y_subset.size() < libsize_t) {
      continue;
    }

    std::vector<double> predictions = SimplexProjection(y_subset, x_subset, nbmat, libsize, E);
    double rho = PearsonCor(y_subset, predictions);
    rhos.push_back(rho);
  }

  return rhos;
}

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
  for (int libsize : libsizes) {
    if (libsize >= E + 1 && libsize <= y_size - 2) {
      libsizevec.push_back(libsize);
    }
  }

  // Sort libsizevec in descending order
  std::sort(libsizevec.begin(), libsizevec.end(), std::greater<int>());

  // Initialize the result vector
  std::vector<std::vector<double>> result;

  // Perform the operations for each libsize in libsizevec using RcppThread
  RcppThread::parallelFor(0, libsizevec.size(), [&](size_t i) {
    int libsize = libsizevec[i];

    std::vector<double> rhos = CalculateRhoForLibsize(y, x, nbmat, libsize, E);

    // Calculate the mean of rhos
    double rhos_mean = CppMean(rhos, true);

    // Store the results in the result vector
    result.push_back({static_cast<double>(libsize), rhos_mean});
  });

  // Calculate significance and confidence interval for each result
  for (size_t i = 0; i < result.size(); ++i) {
    double rho = result[i][1];
    double significance = CppSignificance(rho, y_size);
    std::vector<double> confidence_interval = CppConfidence(rho, y_size);

    result[i].push_back(significance);
    result[i].push_back(confidence_interval[0]);
    result[i].push_back(confidence_interval[1]);
  }

  return result;
}

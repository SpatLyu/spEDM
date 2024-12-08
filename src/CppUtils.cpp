#include <iostream>
#include <vector>
#include <unordered_set>
#include <limits>  // for std::numeric_limits
#include "HelperFuns.h"
#include <Rcpp.h>

// Function to calculate the lagged indices
std::vector<std::vector<int>> CppLaggedIndices(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int lagNum) {
  int n = vec.size();
  std::vector<std::vector<int>> result(n);

  // Handle the case when lagNum is 0
  if (lagNum == 0) {
    for (int i = 0; i < n; ++i) {
      result[i] = {i};
    }
    return result;
  }

  // Handle the case when lagNum is greater than 0
  for (int i = 0; i < n; ++i) {
    std::unordered_set<int> visited;
    std::vector<int> current_neighbors;
    std::vector<int> next_neighbors;

    // Collect 1st level neighbors
    for (int j = 0; j < n; ++j) {
      if (nbmat[i][j] == 1 && i != j) {
        current_neighbors.push_back(j);
        visited.insert(j);
      }
    }

    // Collect neighbors up to lagNum
    for (int l = 1; l < lagNum; ++l) {
      for (int neighbor : current_neighbors) {
        for (int j = 0; j < n; ++j) {
          if (nbmat[neighbor][j] == 1 && i != j && visited.find(j) == visited.end()) {
            next_neighbors.push_back(j);
            visited.insert(j);
          }
        }
      }
      current_neighbors = next_neighbors;
      next_neighbors.clear();
    }

    // Convert set to vector and add to result
    result[i].insert(result[i].end(), visited.begin(), visited.end());

    // If no neighbors found, add NA
    if (result[i].empty()) {
      result[i].push_back(std::numeric_limits<int>::min());
    }
  }

  return result;
}

// Wrapper function to calculate lagged indices and return a List
// [[Rcpp::export]]
Rcpp::List RcppLaggedIndices(const Rcpp::NumericVector& vec,
                             const Rcpp::NumericMatrix& nbmat,
                             int lagNum) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Calculate lagged indices
  std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec_std, nbmat_std, lagNum);

  // Convert std::vector<std::vector<int>> to Rcpp::List
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = Rcpp::wrap(lagged_indices[i]);
  }

  return result;
}

// Function to generate embeddings
std::vector<std::vector<double>> GenEmbeddings(const std::vector<double>& vec,
                                               const std::vector<std::vector<int>>& nbmat,
                                               int E) {
  int n = vec.size();
  std::vector<std::vector<double>> embeddings(n, std::vector<double>(E));

  for (int e = 0; e < E; ++e) {
    int lagNum = E - 1;
    std::vector<std::vector<int>> lagged_indices = CppLaggedIndices(vec, nbmat, lagNum);

    for (int i = 0; i < n; ++i) {
      std::vector<double> lagged_values;
      for (int index : lagged_indices[i]) {
        if (!isNA(index)) {
          lagged_values.push_back(vec[index]);
        }
      }
      embeddings[i][e] = CppMean(lagged_values, true);
    }
  }

  return embeddings;
}

// Wrapper function to generate embeddings and return a NumericMatrix
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppGenEmbeddings(const Rcpp::NumericVector& vec,
                                      const Rcpp::NumericMatrix& nbmat,
                                      int E) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<int>>
  int n = nbmat.nrow();
  std::vector<std::vector<int>> nbmat_std(n, std::vector<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      nbmat_std[i][j] = nbmat(i, j);
    }
  }

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = GenEmbeddings(vec_std, nbmat_std, E);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  int rows = embeddings.size();
  int cols = embeddings[0].size();
  Rcpp::NumericMatrix result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = embeddings[i][j];
    }
  }

  return result;
}

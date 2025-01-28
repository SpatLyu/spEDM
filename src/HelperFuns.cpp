#include <RcppThread.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
unsigned int DetectMaxNumThreads(){
  unsigned int max_threads = std::thread::hardware_concurrency();
  return max_threads;
}

/**
 * Determine the optimal embedding dimension (E) based on the evaluation metrics.
 *
 * This function takes a matrix `Emat` with columns "E", "rho", "mae", and "rmse".
 * It selects the optialmal embedding dimension (E) by first maximizing "rho", then minimizing "rmse",
 * and finally minimizing "mae" if necessary.
 *
 * @param Emat A NumericMatrix with four columns: "E", "rho", "mae", and "rmse".
 * @return The optimal embedding dimension (E) as an integer.
 */
// [[Rcpp::export]]
int OptEmdedDim(Rcpp::NumericMatrix Emat) {
  // Check if the input matrix has exactly 4 columns
  if (Emat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: E, rho, mae, and rmse.");
  }

  // Initialize variables to store the optialmal row index and its metrics
  int optialmal_row = 0;
  double optialmal_rho = Emat(0, 1); // Initialize with the first row's rho
  double optialmal_rmse = Emat(0, 3); // Initialize with the first row's rmse
  double optialmal_mae = Emat(0, 2); // Initialize with the first row's mae

  // Iterate through each row of the matrix
  for (int i = 1; i < Emat.nrow(); ++i) {
    double current_rho = Emat(i, 1); // Current row's rho
    double current_rmse = Emat(i, 3); // Current row's rmse
    double current_mae = Emat(i, 2); // Current row's mae

    // Compare rho values first
    if (current_rho > optialmal_rho) {
      optialmal_row = i;
      optialmal_rho = current_rho;
      optialmal_rmse = current_rmse;
      optialmal_mae = current_mae;
    } else if (current_rho == optialmal_rho) {
      // If rho is equal, compare rmse values
      if (current_rmse < optialmal_rmse) {
        optialmal_row = i;
        optialmal_rho = current_rho;
        optialmal_rmse = current_rmse;
        optialmal_mae = current_mae;
      } else if (current_rmse == optialmal_rmse) {
        // If rmse is also equal, compare mae values
        if (current_mae < optialmal_mae) {
          optialmal_row = i;
          optialmal_rho = current_rho;
          optialmal_rmse = current_rmse;
          optialmal_mae = current_mae;
        }
      }
    }
  }

  // Return the optimal embedding dimension (E) from the optialmal row
  return Emat(optialmal_row, 0);
}

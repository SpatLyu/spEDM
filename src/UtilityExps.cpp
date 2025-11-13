#include <vector>
#include <limits>
#include <cmath>
#include "NumericUtils.h"
#include <RcppThread.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng = false)]]
unsigned int DetectMaxNumThreads(){
  unsigned int max_threads = std::thread::hardware_concurrency();
  return max_threads;
}

/**
 * Determine the optimal embedding dimension (E) and number of nearest neighbors (k).
 *
 * This function selects the best (E, k) combination based on:
 *   1. Maximizing rho
 *   2. Minimizing rmse
 *   3. Minimizing mae
 *   4. If still tied, choosing smallest k, then smallest E
 * A warning is issued when tie-breaking by k and E is used.
 *
 * @param Emat A NumericMatrix with 5 columns: E, k, rho, mae, rmse.
 * @return IntegerVector with optimal E and k.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptEmbedDim(Rcpp::NumericMatrix Emat) {
  if (Emat.ncol() != 5) {
    Rcpp::stop("Input matrix must have exactly 5 columns: E, k, rho, mae, and rmse.");
  }
  if (Emat.nrow() == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  int n = Emat.nrow();
  int opt_row = 0;
  bool used_kE_tiebreak = false;

  struct OptRecord {
    double rho, rmse, mae;
    int k, E;
  };

  OptRecord best{
    Emat(0, 2), // rho
    Emat(0, 4), // rmse
    Emat(0, 3), // mae
    static_cast<int>(Emat(0, 1)), // k
    static_cast<int>(Emat(0, 0))  // E
  };

  for (int i = 1; i < n; ++i) {
    double rho  = Emat(i, 2);
    double mae  = Emat(i, 3);
    double rmse = Emat(i, 4);
    int k = static_cast<int>(Emat(i, 1));
    int E = static_cast<int>(Emat(i, 0));

    bool rho_equal   = doubleNearlyEqual(rho, best.rho);
    bool rmse_equal  = doubleNearlyEqual(rmse, best.rmse);
    bool mae_equal   = doubleNearlyEqual(mae, best.mae);
    // Prevents false positives caused by minimal floating-point deviations
    bool rho_better  = !rho_equal && rho > best.rho;
    bool rmse_better = !rmse_equal && rmse < best.rmse; // smaller is better
    bool mae_better  = !mae_equal && mae < best.mae;    // smaller is better

    if (rho_better ||
        (rho_equal && rmse_better) ||
        (rho_equal && rmse_equal && mae_better)) {
      best = {rho, rmse, mae, k, E};
      opt_row = i;
      used_kE_tiebreak = false;
    } else if (rho_equal && rmse_equal && mae_equal) {
      bool tie_better = (k < best.k) || (k == best.k && E < best.E);
      if (tie_better) {
        best.k = k;
        best.E = E;
        opt_row = i;
        used_kE_tiebreak = true;
      }
    }
  }

  if (used_kE_tiebreak) {
    Rcpp::warning("Tie in evaluation metrics resolved by selecting smallest k, then smallest E.");
  }

  return Rcpp::IntegerVector::create(
    static_cast<int>(Emat(opt_row, 0)), // E
    static_cast<int>(Emat(opt_row, 1))  // k
  );
}

/**
 * Determine the optimal theta parameter based on evaluation metrics.
 *
 * This function takes a NumericMatrix `Thetamat` with columns:
 * "theta", "rho", "mae", and "rmse".
 * The selection criteria are:
 *  - Maximize "rho"
 *  - Minimize "rmse" if "rho" ties
 *  - Minimize "mae" if "rho" and "rmse" tie
 * If multiple rows tie on these metrics (within a tolerance),
 * preference is given to the theta closest to 1.
 * Warnings are issued when tie-breaking occurs or when all metrics are identical.
 *
 * @param Thetamat A NumericMatrix with four columns: theta, rho, mae, and rmse.
 * @return The optimal theta parameter as a double.
 */
// [[Rcpp::export(rng = false)]]
double OptThetaParm(Rcpp::NumericMatrix Thetamat) {
  if (Thetamat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: theta, rho, mae, and rmse.");
  }

  int n = Thetamat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  std::vector<int> best_rows;
  double best_rho  = Thetamat(0, 1);
  double best_mae  = Thetamat(0, 2);
  double best_rmse = Thetamat(0, 3);
  best_rows.push_back(0);

  for (int i = 1; i < n; ++i) {
    double rho  = Thetamat(i, 1);
    double mae  = Thetamat(i, 2);
    double rmse = Thetamat(i, 3);

    bool rho_equal   = doubleNearlyEqual(rho, best_rho);
    bool rmse_equal  = doubleNearlyEqual(rmse, best_rmse);
    bool mae_equal   = doubleNearlyEqual(mae, best_mae);
    // Prevents false positives caused by minimal floating-point deviations
    bool rho_better  = !rho_equal && rho > best_rho;
    bool rmse_better = !rmse_equal && rmse < best_rmse; // smaller is better
    bool mae_better  = !mae_equal && mae < best_mae;    // smaller is better

    if (rho_better ||
        (rho_equal && rmse_better) ||
        (rho_equal && rmse_equal && mae_better)) {
      best_rows.clear();
      best_rows.push_back(i);
      best_rho  = rho;
      best_rmse = rmse;
      best_mae  = mae;
    } else if (rho_equal && rmse_equal && mae_equal) {
      best_rows.push_back(i);
    }
  }

  // If only one best row, return directly
  if (best_rows.size() == 1) {
    return Thetamat(best_rows[0], 0);
  }

  // Check if *all* metrics are identical (within tolerance)
  bool all_equal = true;
  for (size_t i = 1; i < best_rows.size(); ++i) {
    int r0 = best_rows[0];
    int r1 = best_rows[i];
    if (!(doubleNearlyEqual(Thetamat(r0, 1), Thetamat(r1, 1)) &&
          doubleNearlyEqual(Thetamat(r0, 2), Thetamat(r1, 2)) &&
          doubleNearlyEqual(Thetamat(r0, 3), Thetamat(r1, 3)))) {
      all_equal = false;
      break;
    }
  }

  // if *all* metrics are identical, select closest to 1
  double selected_theta = std::numeric_limits<double>::quiet_NaN();
  double min_dist_to_1 = std::numeric_limits<double>::max();

  for (int i : best_rows) {
    double theta = Thetamat(i, 0);
    double dist = std::fabs(theta - 1.0);
    if (dist < min_dist_to_1) {
      min_dist_to_1 = dist;
      selected_theta = theta;
    }
  }

  if (all_equal) {
    Rcpp::warning("All evaluation metrics are identical within tolerance; choosing theta == 1 if available, otherwise closest to 1.");
  } else {
    Rcpp::warning("Tied best evaluation metrics within tolerance; choosing theta == 1 if available, otherwise closest to 1.");
  }

  return selected_theta;
}

/**
 * Select the optimal embedding dimension (E) and number of nearest neighbors (k)
 * from a 4-column matrix: E, k, performance metric, and p-value.
 *
 * Only rows with p-value <= 0.05 are considered.
 * Among them, select the row with:
 *   1. Highest metric (compared using relative tolerance for robustness),
 *   2. If tie, smallest k,
 *   3. If still tie, smallest E.
 *
 * If multiple rows tie on the best metric (within tolerance), a warning is issued
 * and the combination with the smallest k and E is chosen.
 *
 * If no valid rows (p <= 0.05) exist, the function stops with an error.
 *
 * @param Emat NumericMatrix with columns: E, k, metric, and p-value.
 * @return IntegerVector of length 2: optimal E and k.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat) {
  if (Emat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: E, k, metric, and p-value.");
  }

  int n = Emat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  // Filter valid rows with p <= 0.05
  std::vector<int> valid_rows;
  valid_rows.reserve(n);
  for (int i = 0; i < n; ++i) {
    double p = Emat(i, 3);
    if (p < 0.05 || doubleNearlyEqual(p, 0.05)) {
      valid_rows.push_back(i);
    }
  }

  if (valid_rows.empty()) {
    Rcpp::stop("No valid rows with p-value <= 0.05. The chosen neighborhood parameter may be unreasonable or there's no causal relationship. Please consider resetting.");
  }

  int optimal_row = valid_rows[0];
  double best_metric = Emat(optimal_row, 2);
  int best_k = static_cast<int>(Emat(optimal_row, 1));
  int best_E = static_cast<int>(Emat(optimal_row, 0));
  int tie_count = 1;

  for (size_t i = 1; i < valid_rows.size(); ++i) {
    int row = valid_rows[i];
    double current_metric = Emat(row, 2);
    int current_k = static_cast<int>(Emat(row, 1));
    int current_E = static_cast<int>(Emat(row, 0));

    bool metric_equal  = doubleNearlyEqual(current_metric, best_metric);
    bool metric_better = !metric_equal && current_metric > best_metric;

    if (metric_better) {
      optimal_row = row;
      best_metric = current_metric;
      best_k = current_k;
      best_E = current_E;
      tie_count = 1;
    } else if (metric_equal) {
      ++tie_count;
      bool tie_better = (current_k < best_k) ||
      (current_k == best_k && current_E < best_E);
      if (tie_better) {
        optimal_row = row;
        best_k = current_k;
        best_E = current_E;
      }
    }
  }

  if (tie_count > 1) {
    Rcpp::warning("Multiple parameter sets have equal optimal metric; using smallest k and E.");
  }

  return Rcpp::IntegerVector::create(best_E, best_k);
}

/**
 * @title Select Optimal Parameters Based on Causality Metrics
 *
 * @description
 * This function identifies the optimal parameter combination (E, k, tau)
 * from a matrix of evaluation results containing causality measures.
 * The selection is performed hierarchically based on the following criteria:
 *
 * 1. Highest positive causality (`pos`)
 * 2. Highest negative causality (`neg`) when positive causality is equal
 * 3. Highest dark causality (`dark`) when both above are equal
 * 4. In case of complete ties, the smallest E, then smallest tau, then smallest k are chosen.
 *
 * @param Emat A numeric matrix with exactly 6 columns in the order:
 *   E, k, tau, pos, neg, and dark.
 *   Each row represents one evaluation record.
 *
 * @return An integer vector of length 3: the optimal (E, k, tau).
 *
 * @details
 * If multiple parameter sets yield identical causality values (within floating-point tolerance),
 * the function applies a deterministic tie-breaking rule favoring simpler embeddings
 * (smaller E, tau, and k).
 *
 * @note
 * If a tie-break occurs, a warning is issued indicating the decision criteria.
 *
 * @examples
 * \dontrun{
 *   result <- OptPCparms(Emat)
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptPCparms(Rcpp::NumericMatrix Emat) {
  if (Emat.ncol() != 6) {
    Rcpp::stop("Input matrix must have exactly 6 columns: E, k, tau, pos, neg, and dark.");
  }
  if (Emat.nrow() == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  int n = Emat.nrow();
  int opt_row = 0;
  bool used_kE_tiebreak = false;

  struct OptRecord {
    int E, k, tau;
    double pos, neg, dark;
  };

  OptRecord best{
    static_cast<int>(Emat(0, 0)),  // E
    static_cast<int>(Emat(0, 1)),  // k
    static_cast<int>(Emat(0, 2)),  // tau
    Emat(0, 3),                    // pos
    Emat(0, 4),                    // neg
    Emat(0, 5)                     // dark
  };

  for (int i = 1; i < n; ++i) {
    int E = static_cast<int>(Emat(i, 0));
    int k = static_cast<int>(Emat(i, 1));
    int tau = static_cast<int>(Emat(i, 2));
    double pos  = Emat(i, 3);
    double neg  = Emat(i, 4);
    double dark = Emat(i, 5);

    bool pos_equal  = doubleNearlyEqual(pos, best.pos);
    bool neg_equal  = doubleNearlyEqual(neg, best.neg);
    bool dark_equal = doubleNearlyEqual(dark, best.dark);

    bool pos_better  = !pos_equal  && pos  > best.pos;
    bool neg_better  = !neg_equal  && neg  > best.neg;
    bool dark_better = !dark_equal && dark > best.dark;

    if (pos_better ||
        (pos_equal && neg_better) ||
        (pos_equal && neg_equal && dark_better)) {
      best = {E, k, tau, pos, neg, dark};
      opt_row = i;
      used_kE_tiebreak = false;
    } else if (pos_equal && neg_equal && dark_equal) {
      bool tie_better = (E < best.E) ||
        (E == best.E && tau < best.tau) ||
        (E == best.E && tau == best.tau && k < best.k);
      if (tie_better) {
        best.E = E;
        best.tau = tau;
        best.k = k;
        opt_row = i;
        used_kE_tiebreak = true;
      }
    }
  }

  if (used_kE_tiebreak) {
    Rcpp::warning("Tie in evaluation metrics resolved by selecting smallest E, then smallest tau, finally smallest k.");
  }

  return Rcpp::IntegerVector::create(
    static_cast<int>(Emat(opt_row, 0)),  // E
    static_cast<int>(Emat(opt_row, 1)),  // k
    static_cast<int>(Emat(opt_row, 2))   // tau
  );
}

/**
 * This function takes a NumericMatrix as input and returns a matrix
 * containing the row and column indices of all non-NA elements in the input matrix.
 *
 * The processing order can be controlled using the `byrow` parameter:
 *   - If `byrow` is true, the matrix is processed row by row.
 *   - If `byrow` is false, the matrix is processed column by column.
 *
 * Parameters:
 *   - mat: A NumericMatrix object that is to be processed.
 *   - byrow: A boolean parameter to control the processing order.
 *     - If true, the matrix is processed row by row (default is true).
 *     - If false, the matrix is processed column by column.
 *
 * Returns:
 *   - A NumericMatrix with two columns:
 *     - The first column contains the row indices,
 *     - The second column contains the column indices of non-NA elements.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix MatNotNAIndice(Rcpp::NumericMatrix mat, bool byrow = true) {
  // Initialize vectors to store the row and column indices of non-NA elements
  std::vector<double> row_indices;
  std::vector<double> col_indices;

  // Get the number of rows and columns in the input matrix
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  // Loop through the matrix depending on the value of 'byrow'
  if (byrow) {
    // Process by row (row-wise iteration)
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        // Check if the element is not NA
        if (!Rcpp::NumericMatrix::is_na(mat(i, j))) {
          // Record the row and column indices (1-based indexing for R compatibility)
          row_indices.push_back(i + 1);
          col_indices.push_back(j + 1);
        }
      }
    }
  } else {
    // Process by column (column-wise iteration)
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < nrow; i++) {
        // Check if the element is not NA
        if (!Rcpp::NumericMatrix::is_na(mat(i, j))) {
          // Record the row and column indices (1-based indexing for R compatibility)
          row_indices.push_back(i + 1);
          col_indices.push_back(j + 1);
        }
      }
    }
  }

  // Create a NumericMatrix to store the result
  int n = row_indices.size();
  Rcpp::NumericMatrix result(n, 2);

  // Fill the result matrix with the row and column indices
  for (int i = 0; i < n; i++) {
    result(i, 0) = row_indices[i];
    result(i, 1) = col_indices[i];
  }

  // Return the result matrix
  return result;
}

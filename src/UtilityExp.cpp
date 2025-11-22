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
 * Select the optimal embedding parameters (E, k, tau) based on multiple 
 * criteria for simplex projection forecasting.
 *
 * The input matrix must contain six columns in the following order:
 *   1. E      (embedding dimension)
 *   2. k      (number of nearest neighbors)
 *   3. tau    (time or spatial lag)
 *   4. rho    (cross mapping skill, larger is better)
 *   5. mae    (mean absolute error, smaller is better)
 *   6. rmse   (root mean squared error, smaller is better)
 *
 * Selection rules are evaluated in the following order:
 *   1. Maximize rho
 *   2. Minimize rmse
 *   3. Minimize mae
 *   4. If all three metrics are equal within numerical tolerance, resolve ties
 *      by selecting the smallest E, then the smallest tau, then the smallest k.
 *
 * A warning is issued if tie breaking by E, tau, and k is used.
 *
 * @param Emat A NumericMatrix with six columns: E, k, tau, rho, mae, rmse.
 * @return IntegerVector of length three with optimal E, k, and tau in this order.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptSimplexParm(Rcpp::NumericMatrix Emat) {

  if (Emat.ncol() != 6) {
    Rcpp::stop("Input matrix must have exactly six columns: E, k, tau, rho, mae, rmse.");
  }
  if (Emat.nrow() == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  int n = Emat.nrow();
  int opt_row = 0;
  bool used_tiebreak = false;

  struct OptRecord {
    double rho, rmse, mae;
    int E, tau, k;
  };

  OptRecord best{
    Emat(0, 3),                     // rho
    Emat(0, 5),                     // rmse
    Emat(0, 4),                     // mae
    static_cast<int>(Emat(0, 0)),   // E
    static_cast<int>(Emat(0, 2)),   // tau
    static_cast<int>(Emat(0, 1))    // k
  };

  for (int i = 1; i < n; ++i) {

    double rho  = Emat(i, 3);
    double mae  = Emat(i, 4);
    double rmse = Emat(i, 5);

    int E   = static_cast<int>(Emat(i, 0));
    int k   = static_cast<int>(Emat(i, 1));
    int tau = static_cast<int>(Emat(i, 2));

    bool rho_equal  = doubleNearlyEqual(rho, best.rho);
    bool rmse_equal = doubleNearlyEqual(rmse, best.rmse);
    bool mae_equal  = doubleNearlyEqual(mae, best.mae);

    bool rho_better  = !rho_equal && rho > best.rho;
    bool rmse_better = rho_equal && !rmse_equal && rmse < best.rmse;
    bool mae_better  = rho_equal && rmse_equal && !mae_equal && mae < best.mae;

    if (rho_better || rmse_better || mae_better) {
      best = {rho, rmse, mae, E, tau, k};
      opt_row = i;
      used_tiebreak = false;
    }
    else if (rho_equal && rmse_equal && mae_equal) {

      bool tie_better =
        (E < best.E) ||
        (E == best.E && tau < best.tau) ||
        (E == best.E && tau == best.tau && k < best.k);

      if (tie_better) {
        best.E = E;
        best.tau = tau;
        best.k = k;
        opt_row = i;
        used_tiebreak = true;
      }
    }
  }

  if (used_tiebreak) {
    Rcpp::warning("Tie in evaluation metrics was resolved by selecting smallest E, then smallest tau, then smallest k.");
  }

  return Rcpp::IntegerVector::create(
    static_cast<int>(Emat(opt_row, 0)),   // E
    static_cast<int>(Emat(opt_row, 1)),   // k
    static_cast<int>(Emat(opt_row, 2))    // tau
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
 * Select the optimal embedding parameters (E k tau) from intersection cardinality result.
 *
 * The input matrix must contain the following columns in this order:
 *   1. E      embedding dimension
 *   2. k      number of nearest neighbors
 *   3. tau    lag parameter
 *   4. metric performance score to be maximized
 *   5. p      p value used for significance screening
 *
 * Only rows with p value less than or equal to 0.05 are considered valid.
 * Among the valid rows the selection follows these rules:
 *   1. Maximize metric with relative tolerance comparison
 *   2. If metric is equal within tolerance choose smallest E
 *   3. If E is equal choose smallest tau
 *   4. If tau is equal choose smallest k
 *
 * A warning is issued if multiple rows tie on the metric and the final
 * choice is determined by E tau and k.
 *
 * If no valid rows exist the function stops with an error message.
 *
 * @param Emat NumericMatrix with columns: E k tau metric p.
 * @return IntegerVector containing E k tau in this order.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat) {

  if (Emat.ncol() != 5) {
    Rcpp::stop("Input matrix must have exactly five columns: E k tau metric and p value.");
  }

  int n = Emat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  std::vector<int> valid_rows;
  valid_rows.reserve(n);

  for (int i = 0; i < n; ++i) {
    double p = Emat(i, 4);
    if (p < 0.05 || doubleNearlyEqual(p, 0.05)) {
      valid_rows.push_back(i);
    }
  }

  if (valid_rows.empty()) {
    Rcpp::stop("No valid rows with p value less than or equal to 0.05. The chosen neighborhood parameter may be unreasonable or there may be no causal relationship. Consider resetting.");
  }

  int opt_row = valid_rows[0];
  double best_metric = Emat(opt_row, 3);
  int best_E   = static_cast<int>(Emat(opt_row, 0));
  int best_k   = static_cast<int>(Emat(opt_row, 1));
  int best_tau = static_cast<int>(Emat(opt_row, 2));
  int tie_count = 1;

  for (size_t i = 1; i < valid_rows.size(); ++i) {
    int row = valid_rows[i];

    double metric = Emat(row, 3);
    int E   = static_cast<int>(Emat(row, 0));
    int k   = static_cast<int>(Emat(row, 1));
    int tau = static_cast<int>(Emat(row, 2));

    bool metric_equal  = doubleNearlyEqual(metric, best_metric);
    bool metric_better = !metric_equal && metric > best_metric;

    if (metric_better) {
      opt_row = row;
      best_metric = metric;
      best_E = E;
      best_tau = tau;
      best_k = k;
      tie_count = 1;
    }
    else if (metric_equal) {
      tie_count++;

      bool tie_better =
        (E < best_E) ||
        (E == best_E && tau < best_tau) ||
        (E == best_E && tau == best_tau && k < best_k);

      if (tie_better) {
        opt_row = row;
        best_E = E;
        best_tau = tau;
        best_k = k;
      }
    }
  }

  if (tie_count > 1) {
    Rcpp::warning("Multiple parameter sets share the best metric. The final choice was determined by smallest E then smallest tau then smallest k.");
  }

  return Rcpp::IntegerVector::create(best_E, best_k, best_tau);
}

/**
 * @title Select Optimal Parameters Based on Causality Metrics for Pattern Causality
 *
 * @description
 * This function identifies the optimal parameter combination (E, k, tau)
 * from a matrix of evaluation results containing causality measures.
 * Users may choose which causality metric to prioritize using the `maximize`
 * argument. The hierarchical selection rules follow:
 *
 * - maximize = "positive":  pos → dark → neg
 * - maximize = "negative":  neg → dark → pos
 * - maximize = "dark":      dark → pos → neg
 *
 * When all three metrics are tied (within floating-point tolerance),
 * the function applies deterministic tie-breaking based on the smallest
 * E, then tau, then k. A warning is issued when tie-breaking occurs.
 *
 * @param Emat A numeric matrix with 6 columns in the order:
 *   E, k, tau, pos, neg, dark.
 * @param maximize A character string specifying which metric to maximize:
 *   one of "positive", "negative", or "dark".
 *
 * @return An integer vector (E, k, tau) giving the optimal parameters.
 *
 * @examples
 * \dontrun{
 *   OptPCparm(Emat, maximize = "dark")
 * }
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptPCparm(Rcpp::NumericMatrix Emat,
                              std::string maximize = "positive") {

  if (Emat.ncol() != 6) {
    Rcpp::stop("Input matrix must have exactly 6 columns: E, k, tau, pos, neg, dark.");
  }
  if (Emat.nrow() == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  // Validate maximize argument
  if (maximize != "positive" && maximize != "negative" && maximize != "dark") {
    Rcpp::stop("maximize must be one of: 'positive', 'negative', 'dark'.");
  }

  // Determine metric priority order
  // Column indices: pos = 3, neg = 4, dark = 5
  std::vector<int> priority(3);
  if (maximize == "positive") {
    // pos → dark → neg
    priority = {3, 5, 4};
  } else if (maximize == "negative") {
    // neg → dark → pos
    priority = {4, 5, 3};
  } else {
    // maximize == "dark": dark → pos → neg
    priority = {5, 3, 4};
  }

  int n = Emat.nrow();
  int opt_row = 0;
  bool used_tiebreak = false;

  struct OptRecord {
    int E, k, tau;
    double pos, neg, dark;
  };

  // Initialize best record with the first row
  OptRecord best{
    (int)Emat(0, 0),
    (int)Emat(0, 1),
    (int)Emat(0, 2),
    Emat(0, 3),
    Emat(0, 4),
    Emat(0, 5)
  };

  // Helper to map metric index to value
  auto get_metric = [&](const OptRecord &r, int idx) {
    if (idx == 3) return r.pos;
    if (idx == 4) return r.neg;
    return r.dark;  // idx == 5
  };

  for (int i = 1; i < n; ++i) {
    OptRecord cur{
      (int)Emat(i, 0),
      (int)Emat(i, 1),
      (int)Emat(i, 2),
      Emat(i, 3),
      Emat(i, 4),
      Emat(i, 5)
    };

    bool better = false;
    bool equal_all = true;

    // Hierarchical metric comparison
    for (int p = 0; p < 3; ++p) {
      int idx = priority[p];
      double a = get_metric(cur, idx);
      double b = get_metric(best, idx);

      if (!doubleNearlyEqual(a, b)) {
        if (a > b) better = true;
        equal_all = false;
        break;
      }
    }

    if (better) {
      best = cur;
      opt_row = i;
      used_tiebreak = false;
      continue;
    }

    // Complete tie across all three metrics → deterministic tie-breaking
    if (equal_all) {
      bool tie_better =
        (cur.E < best.E) ||
        (cur.E == best.E && cur.tau < best.tau) ||
        (cur.E == best.E && cur.tau == best.tau && cur.k < best.k);

      if (tie_better) {
        best = cur;
        opt_row = i;
        used_tiebreak = true;
      }
    }
  }

  if (used_tiebreak) {
    Rcpp::warning("Tie in evaluation metrics resolved by selecting smallest E, then tau, then k.");
  }

  return Rcpp::IntegerVector::create(
    (int)Emat(opt_row, 0),
    (int)Emat(opt_row, 1),
    (int)Emat(opt_row, 2)
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

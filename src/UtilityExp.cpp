#include <vector>
#include <limits>
#include <cmath>
#include "NumericUtils.h"
#include <thread>
#include <RcppArmadillo.h>

// [[Rcpp::export(rng = false)]]
unsigned int DetectMaxNumThreads(){
  unsigned int max_threads = std::thread::hardware_concurrency();
  return max_threads;
}

/**
 * Select the optimal embedding parameters (E k tau) for simplex projection
 * using a single pass global scan with full tie tracking.
 *
 * The input matrix must contain six columns in the following order:
 *   1. E      embedding dimension
 *   2. k      number of nearest neighbors
 *   3. tau    time or spatial lag
 *   4. rho    cross mapping skill which is maximized
 *   5. mae    mean absolute error which is minimized
 *   6. rmse   root mean squared error which is minimized
 *
 * During the scan the algorithm keeps a list of all globally optimal rows.
 * For each row the evaluation rules are:
 *   1. A row is better if its rho is strictly larger within tolerance
 *   2. If rho is equal a row is better if rmse is smaller within tolerance
 *   3. If rho and rmse are equal a row is better if mae is smaller within tolerance
 *   4. If all three metrics are equal the row is appended to the set of best rows
 *
 * After scanning:
 *   If only one row is optimal its parameters are returned.
 *   If multiple rows are optimal a warning is issued and the final choice is
 *   determined by selecting the smallest E, then the smallest tau, then the smallest k.
 *
 * @param Emat A NumericMatrix with columns: E k tau rho mae rmse.
 * @return IntegerVector of length three that contains E k tau in this order.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptSimplexParm(Rcpp::NumericMatrix Emat) {

  if (Emat.ncol() != 6) {
    Rcpp::stop("Input matrix must have exactly six columns: E k tau rho mae rmse.");
  }
  int n = Emat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  // initialize best record using row 0
  double best_rho  = Emat(0, 3);
  double best_rmse = Emat(0, 5);
  double best_mae  = Emat(0, 4);

  std::vector<int> best_rows;
  best_rows.reserve(n);
  best_rows.push_back(0);

  for (int i = 1; i < n; ++i) {

    double rho  = Emat(i, 3);
    double rmse = Emat(i, 5);
    double mae  = Emat(i, 4);

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
    }
    else if (rho_equal && rmse_equal && mae_equal) {
      best_rows.push_back(i);
    }
  }

  // if only one best row return directly
  if (best_rows.size() == 1) {
    int row = best_rows[0];
    return Rcpp::IntegerVector::create(
      static_cast<int>(Emat(row, 0)),
      static_cast<int>(Emat(row, 1)),
      static_cast<int>(Emat(row, 2))
    );
  }

  // issue tie warning
  Rcpp::warning("Multiple parameter sets share the global optimum for rho rmse and mae. The final choice was determined by smallest E then tau then k.");

  // apply tie breaking rule
  int best_idx = best_rows[0];
  int bestE = Emat(best_idx, 0);
  int bestTau = Emat(best_idx, 2);
  int bestK = Emat(best_idx, 1);

  for (int idx : best_rows) {
    int E = Emat(idx, 0);
    int tau = Emat(idx, 2);
    int k = Emat(idx, 1);

    bool better =
      (E < bestE) ||
      (E == bestE && tau < bestTau) ||
      (E == bestE && tau == bestTau && k < bestK);

    if (better) {
      best_idx = idx;
      bestE = E;
      bestTau = tau;
      bestK = k;
    }
  }

  return Rcpp::IntegerVector::create(bestE, bestK, bestTau);
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
 * Warnings are issued when tie-breaking occurs.
 *
 * @param Thetamat A NumericMatrix with four columns: theta, rho, mae, and rmse.
 * @return The optimal theta parameter as a double.
 */
// [[Rcpp::export(rng = false)]]
double OptThetaParm(Rcpp::NumericMatrix Thetamat) {

  if (Thetamat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly four columns: theta rho mae rmse.");
  }

  int n = Thetamat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  // initialize best metrics using first row
  std::vector<int> best_rows;
  best_rows.push_back(0);

  double best_rho  = Thetamat(0, 1);
  double best_rmse = Thetamat(0, 3);
  double best_mae  = Thetamat(0, 2);

  // global scan through all rows
  for (int i = 1; i < n; ++i) {

    double rho  = Thetamat(i, 1);
    double rmse = Thetamat(i, 3);
    double mae  = Thetamat(i, 2);

    bool rho_equal   = doubleNearlyEqual(rho, best_rho);
    bool rmse_equal  = doubleNearlyEqual(rmse, best_rmse);
    bool mae_equal   = doubleNearlyEqual(mae, best_mae);

    bool rho_better  = (!rho_equal && rho > best_rho);
    bool rmse_better = (rho_equal && !rmse_equal && rmse < best_rmse);
    bool mae_better  = (rho_equal && rmse_equal && !mae_equal && mae < best_mae);

    if (rho_better || rmse_better || mae_better) {
      best_rows.clear();
      best_rows.push_back(i);
      best_rho  = rho;
      best_rmse = rmse;
      best_mae  = mae;
    }
    else if (rho_equal && rmse_equal && mae_equal) {
      best_rows.push_back(i);
    }
  }

  // if only one globally optimal row return directly
  if (best_rows.size() == 1) {
    return Thetamat(best_rows[0], 0);
  }

  // tie exists: choose theta closest to one
  double selected_theta = std::numeric_limits<double>::quiet_NaN();
  double min_dist = std::numeric_limits<double>::max();

  for (int idx : best_rows) {
    double theta = Thetamat(idx, 0);
    double dist = std::fabs(theta - 1.0);

    if (dist < min_dist) {
      min_dist = dist;
      selected_theta = theta;
    }
  }

  Rcpp::warning("Multiple parameter sets share the best evaluation metrics. The final choice is the theta value closest to one.");

  return selected_theta;
}

/**
 * Select the optimal embedding parameters (E k tau) from an intersectional
 * cardinality evaluation matrix using a global scan and full tie tracking.
 *
 * The input matrix must contain the following columns in this order:
 *   1. E      embedding dimension
 *   2. k      number of nearest neighbors
 *   3. tau    lag parameter
 *   4. metric performance score to be maximized
 *   5. p      p value used for significance screening
 *
 * Only rows with p value less than or equal to 0.05 are considered valid.
 *
 * The selection procedure uses a single pass global scan. During the scan
 * the algorithm keeps a list of all rows that jointly achieve the best
 * metric value within numerical tolerance. The comparison rules are:
 *   1. A row is better if its metric is strictly larger within tolerance
 *   2. If the metric is equal a row is appended to the list of best rows
 *
 * After the scan:
 *   If one row is optimal its parameters are returned.
 *   If more than one row is optimal a warning is issued and the final
 *   choice is determined by selecting the smallest E then the smallest tau
 *   then the smallest k.
 *
 * This design prevents misleading warnings that can arise when a local tie
 * appears during the scan but does not correspond to a true global tie.
 *
 * @param Emat A NumericMatrix with columns: E k tau metric p.
 * @return IntegerVector containing E k tau in this order.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat) {

  if (Emat.ncol() != 5) {
    Rcpp::stop("Input matrix must have exactly five columns: E k tau metric p.");
  }
  int n = Emat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  std::vector<int> valid_rows;
  valid_rows.reserve(n);

  // filter rows based on p value
  for (int i = 0; i < n; ++i) {
    double p = Emat(i, 4);
    if (p < 0.05 || doubleNearlyEqual(p, 0.05)) {
      valid_rows.push_back(i);
    }
  }

  if (valid_rows.empty()) {
    Rcpp::stop("No valid rows with p value less than or equal to 0.05. The neighborhood parameter may be unreasonable or there may be no causal relationship.");
  }

  // initialize best set using the first valid row
  int first = valid_rows[0];
  double best_metric = Emat(first, 3);

  std::vector<int> best_rows;
  best_rows.push_back(first);

  // global scan across valid rows
  for (size_t i = 1; i < valid_rows.size(); ++i) {
    int row = valid_rows[i];
    double metric = Emat(row, 3);

    bool metric_equal  = doubleNearlyEqual(metric, best_metric);
    bool metric_better = (!metric_equal && metric > best_metric);

    if (metric_better) {
      best_rows.clear();
      best_rows.push_back(row);
      best_metric = metric;
    }
    else if (metric_equal) {
      best_rows.push_back(row);
    }
  }

  // if only one row is globally optimal return directly
  if (best_rows.size() == 1) {
    int row = best_rows[0];
    return Rcpp::IntegerVector::create(
      static_cast<int>(Emat(row, 0)),
      static_cast<int>(Emat(row, 1)),
      static_cast<int>(Emat(row, 2))
    );
  }

  // issue tie warning
  Rcpp::warning("Multiple parameter sets share the best metric. The final choice is determined by smallest E then smallest tau then smallest k.");

  // tie breaking: smallest E then tau then k
  int best_idx = best_rows[0];
  int bestE = Emat(best_idx, 0);
  int bestTau = Emat(best_idx, 2);
  int bestK = Emat(best_idx, 1);

  for (int idx : best_rows) {
    int E = Emat(idx, 0);
    int tau = Emat(idx, 2);
    int k = Emat(idx, 1);

    bool better =
      (E < bestE) ||
      (E == bestE && tau < bestTau) ||
      (E == bestE && tau == bestTau && k < bestK);

    if (better) {
      best_idx = idx;
      bestE = E;
      bestTau = tau;
      bestK = k;
    }
  }

  return Rcpp::IntegerVector::create(bestE, bestK, bestTau);
}

/**
 * Select optimal embedding parameters (E k tau) for pattern causality metrics
 * using a global scan and full tie tracking across three causality measures.
 *
 * The input matrix must contain six columns in the following order:
 *   1. E      embedding dimension
 *   2. k      number of nearest neighbors
 *   3. tau    lag parameter
 *   4. pos    positive causality score
 *   5. neg    negative causality score
 *   6. dark   dark causality score
 *
 * The argument `maximize` specifies which causality metric should be given
 * primary priority. The three metrics are compared in a hierarchical order:
 *
 *   maximize = "positive":  pos then dark then neg
 *   maximize = "negative":  neg then dark then pos
 *   maximize = "dark":      dark then pos then neg
 *
 * The function performs a single pass global scan. During the scan it collects
 * all rows that achieve the joint optimum across the three prioritized metrics
 * within numerical tolerance. When multiple rows tie for the optimum a final
 * deterministic choice is made by selecting the smallest E then the smallest
 * tau then the smallest k. A warning is issued when this final tie breaking
 * procedure is required.
 *
 * @param Emat A numeric matrix with columns: E k tau pos neg dark.
 * @param maximize A string specifying which metric to prioritize.
 * @return IntegerVector containing E k tau in this order.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptPCparm(Rcpp::NumericMatrix Emat,
                              std::string maximize = "positive") {

  if (Emat.ncol() != 6) {
    Rcpp::stop("Input matrix must have exactly six columns: E k tau pos neg dark.");
  }
  int n = Emat.nrow();
  if (n == 0) {
    Rcpp::stop("Input matrix must not be empty.");
  }

  if (maximize != "positive" && maximize != "negative" && maximize != "dark") {
    Rcpp::stop("maximize must be one of positive negative or dark.");
  }

  // establish metric priority order
  std::vector<int> priority(3);
  if (maximize == "positive") {
    priority = {3, 5, 4};  // pos dark neg
  } else if (maximize == "negative") {
    priority = {4, 5, 3};  // neg dark pos
  } else {
    priority = {5, 3, 4};  // dark pos neg
  }

  // helper to get metric value
  auto get_metric = [&](int row, int idx) {
    return Emat(row, idx);
  };

  // record of the best metrics for comparison
  std::vector<double> best_vals(3);
  for (int j = 0; j < 3; ++j) {
    best_vals[j] = get_metric(0, priority[j]);
  }

  std::vector<int> best_rows;
  best_rows.push_back(0);

  // global scan
  for (int i = 1; i < n; ++i) {

    bool better = false;
    bool equal_all = true;

    for (int p = 0; p < 3; ++p) {
      double a = get_metric(i, priority[p]);
      double b = best_vals[p];

      if (!doubleNearlyEqual(a, b)) {
        if (a > b) {
          better = true;
        }
        equal_all = false;
        break;
      }
    }

    if (better) {
      best_rows.clear();
      best_rows.push_back(i);
      for (int p = 0; p < 3; ++p) {
        best_vals[p] = get_metric(i, priority[p]);
      }
    }
    else if (equal_all) {
      best_rows.push_back(i);
    }
  }

  // if only one globally optimal row return directly
  if (best_rows.size() == 1) {
    int row = best_rows[0];
    return Rcpp::IntegerVector::create(
      static_cast<int>(Emat(row, 0)),
      static_cast<int>(Emat(row, 1)),
      static_cast<int>(Emat(row, 2))
    );
  }

  // issue tie warning
  Rcpp::warning("Multiple parameter sets share the global optimum. The final choice is determined by smallest E then tau then k.");

  // tie break by E then tau then k
  int best_idx = best_rows[0];
  int bestE = Emat(best_idx, 0);
  int bestTau = Emat(best_idx, 2);
  int bestK = Emat(best_idx, 1);

  for (int idx : best_rows) {
    int E = Emat(idx, 0);
    int tau = Emat(idx, 2);
    int k = Emat(idx, 1);

    bool better =
      (E < bestE) ||
      (E == bestE && tau < bestTau) ||
      (E == bestE && tau == bestTau && k < bestK);

    if (better) {
      best_idx = idx;
      bestE = E;
      bestTau = tau;
      bestK = k;
    }
  }

  return Rcpp::IntegerVector::create(bestE, bestK, bestTau);
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

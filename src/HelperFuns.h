#ifndef HelperFuns_H
#define HelperFuns_H

#include <vector>
#include <RcppThread.h>
#include <RcppArmadillo.h>

/**
 * Determine the optimal embedding dimension (E) and number of nearest neighbors (k).
 *
 * This function takes a matrix `Emat` with columns "E", "k", "rho", "mae", and "rmse".
 * It selects the optialmal embedding dimension (E) and number of nearest neighbors (k)
 * by first maximizing "rho", then minimizing "rmse", and finally minimizing "mae" if necessary.
 *
 * @param Emat A NumericMatrix with five columns: "E", k", "rho", "mae", and "rmse".
 * @return The optimal embedding dimension (E) and number of nearest neighbors (k) as an integer vector.
 */
Rcpp::IntegerVector OptEmbedDim(Rcpp::NumericMatrix Emat);

/**
 * Determine the optimal theta parameter based on the evaluation metrics.
 *
 * This function takes a matrix `Thetamat` with columns "theta", "rho", "mae", and "rmse".
 * It selects the optimal theta parameter by first maximizing "rho",
 * then minimizing "rmse", and finally minimizing "mae" if necessary.
 *
 * @param Thetamat A NumericMatrix with four columns: "theta", "rho", "mae", and "rmse".
 * @return The optimal theta parameter as a double.
 */
double OptThetaParm(Rcpp::NumericMatrix Thetamat);

/**
 * Select the optimal embedding dimension (E) and number of nearest neighbors (k)
 * from a 4-column matrix: E, k, performance metric, and p-value.
 * Only rows with p-value <= 0.05 are considered. Among them, rows are compared by:
 *  - highest metric,
 *  - smallest k (tie-breaker),
 *  - smallest E (final tie-breaker).
 *
 * @param Emat A NumericMatrix with four columns: E, k, metric, and p-value.
 * @return IntegerVector of length 2: the optimal E and k.
 */
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat);

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
Rcpp::NumericMatrix MatNotNAIndice(Rcpp::NumericMatrix mat, bool byrow = true);

#endif // HelperFuns

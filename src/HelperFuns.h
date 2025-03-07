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
 * This function takes a NumericMatrix as input and returns a matrix
 * containing the row and column indices of all non-NA elements in the input matrix.
 *
 * Parameters:
 *   - mat: A NumericMatrix object that is to be processed.
 *
 * Returns:
 *   - A NumericMatrix with two columns:
 *     - The first column contains the row indices,
 *     - The second column contains the column indices of non-NA elements.
 */
Rcpp::NumericMatrix MatNotNAIndice(Rcpp::NumericMatrix mat)

#endif // HelperFuns

#ifndef HelperFuns_H
#define HelperFuns_H

#include <RcppThread.h>
#include <RcppArmadillo.h>

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
int OptEmdedDim(Rcpp::NumericMatrix Emat);

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

#endif // HelperFuns

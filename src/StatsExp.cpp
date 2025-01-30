#include <vector>
#include "CppStats.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double RcppMean(const Rcpp::NumericVector& vec,
                bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppMean(y, NA_rm);
}

// [[Rcpp::export]]
double RcppSum(const Rcpp::NumericVector& vec,
               bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppSum(y, NA_rm);
}

// [[Rcpp::export]]
double RcppVariance(const Rcpp::NumericVector& vec,
                    bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);

  // Call the ArmaPearsonCor function
  return CppVariance(y, NA_rm);
}

// [[Rcpp::export]]
double RcppCovariance(const Rcpp::NumericVector& vec1,
                      const Rcpp::NumericVector& vec2,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppMAE function
  return CppCovariance(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
double RcppMAE(const Rcpp::NumericVector& vec1,
               const Rcpp::NumericVector& vec2,
               bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppMAE function
  return CppMAE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
double RcppRMSE(const Rcpp::NumericVector& vec1,
                const Rcpp::NumericVector& vec2,
                bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppRMSE function
  return CppRMSE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export]]
Rcpp::NumericVector RcppAbs(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec1_cpp = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> vec2_cpp = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppAbs function
  std::vector<double> result = CppAbs(vec1_cpp, vec2_cpp);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
Rcpp::NumericVector RcppSumNormalize(const Rcpp::NumericVector& vec, bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_cpp = Rcpp::as<std::vector<double>>(vec);

  // Call the CppSumNormalize function
  std::vector<double> result = CppSumNormalize(vec_cpp, NA_rm);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
double RcppPearsonCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_vec = Rcpp::as<std::vector<double>>(y);
  std::vector<double> y_hat_vec = Rcpp::as<std::vector<double>>(y_hat);

  // Call the ArmaPearsonCor function
  return PearsonCor(y_vec, y_hat_vec, NA_rm);
}

// Rcpp wrapper for PartialCor function
// [[Rcpp::export]]
double RcppPartialCor(Rcpp::NumericVector y,
                      Rcpp::NumericVector y_hat,
                      Rcpp::NumericMatrix controls,
                      bool NA_rm = false,
                      bool linear = false) {

  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> std_controls(controls.ncol());
  for (int i = 0; i < controls.ncol(); ++i) {
    Rcpp::NumericVector covvar = controls.column(i);
    std_controls[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Call the PartialCor function
  return PartialCor(std_y, std_y_hat, std_controls, NA_rm, linear);
}

// [[Rcpp::export]]
double RcppPartialCorTrivar(Rcpp::NumericVector y,
                            Rcpp::NumericVector y_hat,
                            Rcpp::NumericVector control,
                            bool NA_rm = false,
                            bool linear = false) {

  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);
  std::vector<double> std_control = Rcpp::as<std::vector<double>>(control);

  // Call the PartialCor function
  return PartialCorTrivar(std_y, std_y_hat, std_control, NA_rm, linear);
}

// Wrapper function to calculate the confidence interval for a correlation coefficient and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppConfidence(double r, int n, double level = 0.05) {
  // Calculate the confidence interval
  std::vector<double> result = CppConfidence(r, n, level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to perform Linear Trend Removal and return a NumericVector
// [[Rcpp::export]]
Rcpp::NumericVector RcppLinearTrendRM(const Rcpp::NumericVector& vec,
                                      const Rcpp::NumericVector& xcoord,
                                      const Rcpp::NumericVector& ycoord,
                                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
  std::vector<double> xcoord_std = Rcpp::as<std::vector<double>>(xcoord);
  std::vector<double> ycoord_std = Rcpp::as<std::vector<double>>(ycoord);

  // Perform Linear Trend Removal
  std::vector<double> result = LinearTrendRM(vec_std, xcoord_std, ycoord_std, NA_rm);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Rcpp wrapper function for CppSVD
// [[Rcpp::export]]
Rcpp::List RcppSVD(const Rcpp::NumericMatrix& X) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t m = X.nrow();
  size_t n = X.ncol();
  std::vector<std::vector<double>> X_vec(m, std::vector<double>(n));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      X_vec[i][j] = X(i, j);
    }
  }

  // Call the original CppSVD function
  std::vector<std::vector<std::vector<double>>> result = CppSVD(X_vec);

  // Extract results from CppSVD output
  std::vector<std::vector<double>> u = result[0]; // Left singular vectors
  std::vector<double> d = result[1][0];           // Singular values
  std::vector<std::vector<double>> v = result[2]; // Right singular vectors

  // Convert std::vector results to Rcpp objects
  Rcpp::NumericMatrix u_rcpp(m, m);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      u_rcpp(i, j) = u[i][j];
    }
  }

  Rcpp::NumericVector d_rcpp(d.size());
  for (size_t i = 0; i < d.size(); ++i) {
    d_rcpp(i) = d[i];
  }

  Rcpp::NumericMatrix v_rcpp(v.size(), v[0].size());
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v_rcpp(i, j) = v[i][j];
    }
  }

  // Return results as an Rcpp::List to match R's svd() output
  return Rcpp::List::create(
    Rcpp::Named("u") = u_rcpp, // Left singular vectors
    Rcpp::Named("d") = d_rcpp, // Singular values
    Rcpp::Named("v") = v_rcpp  // Right singular vectors
  );
}

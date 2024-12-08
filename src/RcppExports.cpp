// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppLaggedIndices
Rcpp::List RcppLaggedIndices(const Rcpp::NumericVector& vec, const Rcpp::NumericMatrix& nbmat, int lagNum);
RcppExport SEXP _spEDM_RcppLaggedIndices(SEXP vecSEXP, SEXP nbmatSEXP, SEXP lagNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type nbmat(nbmatSEXP);
    Rcpp::traits::input_parameter< int >::type lagNum(lagNumSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLaggedIndices(vec, nbmat, lagNum));
    return rcpp_result_gen;
END_RCPP
}
// RcppGenEmbeddings
Rcpp::NumericMatrix RcppGenEmbeddings(const Rcpp::NumericVector& vec, const Rcpp::NumericMatrix& nbmat, int E);
RcppExport SEXP _spEDM_RcppGenEmbeddings(SEXP vecSEXP, SEXP nbmatSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type nbmat(nbmatSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGenEmbeddings(vec, nbmat, E));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _spEDM_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spEDM_RcppLaggedIndices", (DL_FUNC) &_spEDM_RcppLaggedIndices, 3},
    {"_spEDM_RcppGenEmbeddings", (DL_FUNC) &_spEDM_RcppGenEmbeddings, 3},
    {"_spEDM_timesTwo", (DL_FUNC) &_spEDM_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_spEDM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

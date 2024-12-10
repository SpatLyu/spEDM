// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppConfidence
Rcpp::NumericVector RcppConfidence(double r, int n, double level);
RcppExport SEXP _spEDM_RcppConfidence(SEXP rSEXP, SEXP nSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppConfidence(r, n, level));
    return rcpp_result_gen;
END_RCPP
}
// RcppLinearTrendRM
Rcpp::NumericVector RcppLinearTrendRM(const Rcpp::NumericVector& vec, const Rcpp::NumericVector& xcoord, const Rcpp::NumericVector& ycoord, bool NA_rm);
RcppExport SEXP _spEDM_RcppLinearTrendRM(SEXP vecSEXP, SEXP xcoordSEXP, SEXP ycoordSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xcoord(xcoordSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ycoord(ycoordSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLinearTrendRM(vec, xcoord, ycoord, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
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
// RcppGCCMLattice
Rcpp::NumericMatrix RcppGCCMLattice(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& nbmat, const Rcpp::IntegerVector& libsizes, int E);
RcppExport SEXP _spEDM_RcppGCCMLattice(SEXP xSEXP, SEXP ySEXP, SEXP nbmatSEXP, SEXP libsizesSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type nbmat(nbmatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type libsizes(libsizesSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGCCMLattice(x, y, nbmat, libsizes, E));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spEDM_RcppConfidence", (DL_FUNC) &_spEDM_RcppConfidence, 3},
    {"_spEDM_RcppLinearTrendRM", (DL_FUNC) &_spEDM_RcppLinearTrendRM, 4},
    {"_spEDM_RcppLaggedIndices", (DL_FUNC) &_spEDM_RcppLaggedIndices, 3},
    {"_spEDM_RcppGenEmbeddings", (DL_FUNC) &_spEDM_RcppGenEmbeddings, 3},
    {"_spEDM_RcppGCCMLattice", (DL_FUNC) &_spEDM_RcppGCCMLattice, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spEDM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

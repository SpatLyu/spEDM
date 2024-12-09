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
// RcppDist
Rcpp::NumericMatrix RcppDist(const Rcpp::NumericMatrix& matrix);
RcppExport SEXP _spEDM_RcppDist(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppDist(matrix));
    return rcpp_result_gen;
END_RCPP
}
// RcppClosestIndices
Rcpp::IntegerMatrix RcppClosestIndices(const Rcpp::NumericMatrix& distmat, int libsize);
RcppExport SEXP _spEDM_RcppClosestIndices(SEXP distmatSEXP, SEXP libsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< int >::type libsize(libsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppClosestIndices(distmat, libsize));
    return rcpp_result_gen;
END_RCPP
}
// RcppCCMWeight
Rcpp::NumericMatrix RcppCCMWeight(const Rcpp::NumericMatrix& distmat, const Rcpp::IntegerMatrix& closestIndices, int libsize);
RcppExport SEXP _spEDM_RcppCCMWeight(SEXP distmatSEXP, SEXP closestIndicesSEXP, SEXP libsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type closestIndices(closestIndicesSEXP);
    Rcpp::traits::input_parameter< int >::type libsize(libsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppCCMWeight(distmat, closestIndices, libsize));
    return rcpp_result_gen;
END_RCPP
}
// RcppSimplexProjection
Rcpp::NumericVector RcppSimplexProjection(const Rcpp::NumericVector& y, const Rcpp::NumericVector& x, const Rcpp::NumericMatrix& nbmat, int libsize, int E);
RcppExport SEXP _spEDM_RcppSimplexProjection(SEXP ySEXP, SEXP xSEXP, SEXP nbmatSEXP, SEXP libsizeSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type nbmat(nbmatSEXP);
    Rcpp::traits::input_parameter< int >::type libsize(libsizeSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSimplexProjection(y, x, nbmat, libsize, E));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spEDM_RcppConfidence", (DL_FUNC) &_spEDM_RcppConfidence, 3},
    {"_spEDM_RcppLaggedIndices", (DL_FUNC) &_spEDM_RcppLaggedIndices, 3},
    {"_spEDM_RcppGenEmbeddings", (DL_FUNC) &_spEDM_RcppGenEmbeddings, 3},
    {"_spEDM_RcppDist", (DL_FUNC) &_spEDM_RcppDist, 1},
    {"_spEDM_RcppClosestIndices", (DL_FUNC) &_spEDM_RcppClosestIndices, 2},
    {"_spEDM_RcppCCMWeight", (DL_FUNC) &_spEDM_RcppCCMWeight, 3},
    {"_spEDM_RcppSimplexProjection", (DL_FUNC) &_spEDM_RcppSimplexProjection, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spEDM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

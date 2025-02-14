// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppSimplexForecast
Rcpp::NumericVector RcppSimplexForecast(const Rcpp::NumericMatrix& embedding, const Rcpp::NumericVector& target, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, const int& num_neighbors);
RcppExport SEXP _spEDM_RcppSimplexForecast(SEXP embeddingSEXP, SEXP targetSEXP, SEXP libSEXP, SEXP predSEXP, SEXP num_neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type embedding(embeddingSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_neighbors(num_neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSimplexForecast(embedding, target, lib, pred, num_neighbors));
    return rcpp_result_gen;
END_RCPP
}
// RcppSMapForecast
Rcpp::NumericVector RcppSMapForecast(const Rcpp::NumericMatrix& embedding, const Rcpp::NumericVector& target, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, const int& num_neighbors, const double& theta);
RcppExport SEXP _spEDM_RcppSMapForecast(SEXP embeddingSEXP, SEXP targetSEXP, SEXP libSEXP, SEXP predSEXP, SEXP num_neighborsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type embedding(embeddingSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_neighbors(num_neighborsSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSMapForecast(embedding, target, lib, pred, num_neighbors, theta));
    return rcpp_result_gen;
END_RCPP
}
// RcppLaggedVar4Grid
Rcpp::NumericMatrix RcppLaggedVar4Grid(const Rcpp::NumericMatrix& mat, int lagNum);
RcppExport SEXP _spEDM_RcppLaggedVar4Grid(SEXP matSEXP, SEXP lagNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type lagNum(lagNumSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLaggedVar4Grid(mat, lagNum));
    return rcpp_result_gen;
END_RCPP
}
// RcppGenGridEmbeddings
Rcpp::NumericMatrix RcppGenGridEmbeddings(const Rcpp::NumericMatrix& mat, int E, int tau);
RcppExport SEXP _spEDM_RcppGenGridEmbeddings(SEXP matSEXP, SEXP ESEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGenGridEmbeddings(mat, E, tau));
    return rcpp_result_gen;
END_RCPP
}
// RcppLocateGridIndices
int RcppLocateGridIndices(int curRow, int curCol, int totalRow, int totalCol);
RcppExport SEXP _spEDM_RcppLocateGridIndices(SEXP curRowSEXP, SEXP curColSEXP, SEXP totalRowSEXP, SEXP totalColSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type curRow(curRowSEXP);
    Rcpp::traits::input_parameter< int >::type curCol(curColSEXP);
    Rcpp::traits::input_parameter< int >::type totalRow(totalRowSEXP);
    Rcpp::traits::input_parameter< int >::type totalCol(totalColSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLocateGridIndices(curRow, curCol, totalRow, totalCol));
    return rcpp_result_gen;
END_RCPP
}
// RcppSimplex4Grid
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& mat, const Rcpp::IntegerMatrix& lib, const Rcpp::IntegerMatrix& pred, const Rcpp::IntegerVector& E, int tau, int b, int threads);
RcppExport SEXP _spEDM_RcppSimplex4Grid(SEXP matSEXP, SEXP libSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSimplex4Grid(mat, lib, pred, E, tau, b, threads));
    return rcpp_result_gen;
END_RCPP
}
// RcppSMap4Grid
Rcpp::NumericMatrix RcppSMap4Grid(const Rcpp::NumericMatrix& mat, const Rcpp::IntegerMatrix& lib, const Rcpp::IntegerMatrix& pred, const Rcpp::NumericVector& theta, int E, int tau, int b, int threads);
RcppExport SEXP _spEDM_RcppSMap4Grid(SEXP matSEXP, SEXP libSEXP, SEXP predSEXP, SEXP thetaSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSMap4Grid(mat, lib, pred, theta, E, tau, b, threads));
    return rcpp_result_gen;
END_RCPP
}
// RcppGCCM4Grid
Rcpp::NumericMatrix RcppGCCM4Grid(const Rcpp::NumericMatrix& xMatrix, const Rcpp::NumericMatrix& yMatrix, const Rcpp::IntegerVector& libsizes, const Rcpp::IntegerMatrix& pred, int E, int tau, int b, bool simplex, double theta, int threads, bool progressbar);
RcppExport SEXP _spEDM_RcppGCCM4Grid(SEXP xMatrixSEXP, SEXP yMatrixSEXP, SEXP libsizesSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP simplexSEXP, SEXP thetaSEXP, SEXP threadsSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type xMatrix(xMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type yMatrix(yMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type libsizes(libsizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type simplex(simplexSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGCCM4Grid(xMatrix, yMatrix, libsizes, pred, E, tau, b, simplex, theta, threads, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// RcppSCPCM4Grid
Rcpp::NumericMatrix RcppSCPCM4Grid(const Rcpp::NumericMatrix& xMatrix, const Rcpp::NumericMatrix& yMatrix, const Rcpp::NumericMatrix& zMatrix, const Rcpp::IntegerVector& libsizes, const Rcpp::IntegerVector& E, const Rcpp::IntegerVector& tau, const Rcpp::IntegerMatrix& pred, int b, bool simplex, double theta, int threads, bool cumulate, bool progressbar);
RcppExport SEXP _spEDM_RcppSCPCM4Grid(SEXP xMatrixSEXP, SEXP yMatrixSEXP, SEXP zMatrixSEXP, SEXP libsizesSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP predSEXP, SEXP bSEXP, SEXP simplexSEXP, SEXP thetaSEXP, SEXP threadsSEXP, SEXP cumulateSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type xMatrix(xMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type yMatrix(yMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type zMatrix(zMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type libsizes(libsizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type simplex(simplexSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type cumulate(cumulateSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSCPCM4Grid(xMatrix, yMatrix, zMatrix, libsizes, E, tau, pred, b, simplex, theta, threads, cumulate, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// RcppGCMC4Grid
Rcpp::NumericVector RcppGCMC4Grid(const Rcpp::NumericMatrix& xMatrix, const Rcpp::NumericMatrix& yMatrix, const Rcpp::IntegerMatrix& pred, const Rcpp::IntegerVector& E, const Rcpp::IntegerVector& tau, const Rcpp::IntegerVector& b, const Rcpp::IntegerVector& max_r, int threads, bool progressbar);
RcppExport SEXP _spEDM_RcppGCMC4Grid(SEXP xMatrixSEXP, SEXP yMatrixSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP max_rSEXP, SEXP threadsSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type xMatrix(xMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type yMatrix(yMatrixSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type max_r(max_rSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGCMC4Grid(xMatrix, yMatrix, pred, E, tau, b, max_r, threads, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// DetectMaxNumThreads
unsigned int DetectMaxNumThreads();
RcppExport SEXP _spEDM_DetectMaxNumThreads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(DetectMaxNumThreads());
    return rcpp_result_gen;
END_RCPP
}
// OptEmdedDim
int OptEmdedDim(Rcpp::NumericMatrix Emat);
RcppExport SEXP _spEDM_OptEmdedDim(SEXP EmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Emat(EmatSEXP);
    rcpp_result_gen = Rcpp::wrap(OptEmdedDim(Emat));
    return rcpp_result_gen;
END_RCPP
}
// OptThetaParm
double OptThetaParm(Rcpp::NumericMatrix Thetamat);
RcppExport SEXP _spEDM_OptThetaParm(SEXP ThetamatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Thetamat(ThetamatSEXP);
    rcpp_result_gen = Rcpp::wrap(OptThetaParm(Thetamat));
    return rcpp_result_gen;
END_RCPP
}
// RcppLaggedNeighbor4Lattice
Rcpp::List RcppLaggedNeighbor4Lattice(const Rcpp::List& nb, int lagNum);
RcppExport SEXP _spEDM_RcppLaggedNeighbor4Lattice(SEXP nbSEXP, SEXP lagNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< int >::type lagNum(lagNumSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLaggedNeighbor4Lattice(nb, lagNum));
    return rcpp_result_gen;
END_RCPP
}
// RcppLaggedVar4Lattice
Rcpp::List RcppLaggedVar4Lattice(const Rcpp::List& nb, int lagNum);
RcppExport SEXP _spEDM_RcppLaggedVar4Lattice(SEXP nbSEXP, SEXP lagNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< int >::type lagNum(lagNumSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppLaggedVar4Lattice(nb, lagNum));
    return rcpp_result_gen;
END_RCPP
}
// RcppGenLatticeEmbeddings
Rcpp::NumericMatrix RcppGenLatticeEmbeddings(const Rcpp::NumericVector& vec, const Rcpp::List& nb, int E, int tau);
RcppExport SEXP _spEDM_RcppGenLatticeEmbeddings(SEXP vecSEXP, SEXP nbSEXP, SEXP ESEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGenLatticeEmbeddings(vec, nb, E, tau));
    return rcpp_result_gen;
END_RCPP
}
// RcppSimplex4Lattice
Rcpp::NumericMatrix RcppSimplex4Lattice(const Rcpp::NumericVector& x, const Rcpp::List& nb, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, const Rcpp::IntegerVector& E, int tau, int b, int threads);
RcppExport SEXP _spEDM_RcppSimplex4Lattice(SEXP xSEXP, SEXP nbSEXP, SEXP libSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSimplex4Lattice(x, nb, lib, pred, E, tau, b, threads));
    return rcpp_result_gen;
END_RCPP
}
// RcppSMap4Lattice
Rcpp::NumericMatrix RcppSMap4Lattice(const Rcpp::NumericVector& x, const Rcpp::List& nb, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, const Rcpp::NumericVector& theta, int E, int tau, int b, int threads);
RcppExport SEXP _spEDM_RcppSMap4Lattice(SEXP xSEXP, SEXP nbSEXP, SEXP libSEXP, SEXP predSEXP, SEXP thetaSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSMap4Lattice(x, nb, lib, pred, theta, E, tau, b, threads));
    return rcpp_result_gen;
END_RCPP
}
// RcppGCCM4Lattice
Rcpp::NumericMatrix RcppGCCM4Lattice(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& nb, const Rcpp::IntegerVector& libsizes, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, int E, int tau, int b, bool simplex, double theta, int threads, bool progressbar);
RcppExport SEXP _spEDM_RcppGCCM4Lattice(SEXP xSEXP, SEXP ySEXP, SEXP nbSEXP, SEXP libsizesSEXP, SEXP libSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP simplexSEXP, SEXP thetaSEXP, SEXP threadsSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type libsizes(libsizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type simplex(simplexSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGCCM4Lattice(x, y, nb, libsizes, lib, pred, E, tau, b, simplex, theta, threads, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// RcppSCPCM4Lattice
Rcpp::NumericMatrix RcppSCPCM4Lattice(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& z, const Rcpp::List& nb, const Rcpp::IntegerVector& libsizes, const Rcpp::IntegerVector& lib, const Rcpp::IntegerVector& pred, const Rcpp::IntegerVector& E, const Rcpp::IntegerVector& tau, int b, bool simplex, double theta, int threads, bool cumulate, bool progressbar);
RcppExport SEXP _spEDM_RcppSCPCM4Lattice(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP nbSEXP, SEXP libsizesSEXP, SEXP libSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP simplexSEXP, SEXP thetaSEXP, SEXP threadsSEXP, SEXP cumulateSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type libsizes(libsizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lib(libSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type simplex(simplexSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type cumulate(cumulateSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSCPCM4Lattice(x, y, z, nb, libsizes, lib, pred, E, tau, b, simplex, theta, threads, cumulate, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// RcppGCMC4Lattice
Rcpp::NumericVector RcppGCMC4Lattice(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y, const Rcpp::List& nb, const Rcpp::IntegerVector& pred, const Rcpp::IntegerVector& E, const Rcpp::IntegerVector& tau, const Rcpp::IntegerVector& b, const Rcpp::IntegerVector& max_r, int threads, bool progressbar);
RcppExport SEXP _spEDM_RcppGCMC4Lattice(SEXP xSEXP, SEXP ySEXP, SEXP nbSEXP, SEXP predSEXP, SEXP ESEXP, SEXP tauSEXP, SEXP bSEXP, SEXP max_rSEXP, SEXP threadsSEXP, SEXP progressbarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type E(ESEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type max_r(max_rSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progressbar(progressbarSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppGCMC4Lattice(x, y, nb, pred, E, tau, b, max_r, threads, progressbar));
    return rcpp_result_gen;
END_RCPP
}
// RcppMean
double RcppMean(const Rcpp::NumericVector& vec, bool NA_rm);
RcppExport SEXP _spEDM_RcppMean(SEXP vecSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppMean(vec, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppSum
double RcppSum(const Rcpp::NumericVector& vec, bool NA_rm);
RcppExport SEXP _spEDM_RcppSum(SEXP vecSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSum(vec, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppVariance
double RcppVariance(const Rcpp::NumericVector& vec, bool NA_rm);
RcppExport SEXP _spEDM_RcppVariance(SEXP vecSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppVariance(vec, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppCovariance
double RcppCovariance(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2, bool NA_rm);
RcppExport SEXP _spEDM_RcppCovariance(SEXP vec1SEXP, SEXP vec2SEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppCovariance(vec1, vec2, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppMAE
double RcppMAE(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2, bool NA_rm);
RcppExport SEXP _spEDM_RcppMAE(SEXP vec1SEXP, SEXP vec2SEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppMAE(vec1, vec2, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppRMSE
double RcppRMSE(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2, bool NA_rm);
RcppExport SEXP _spEDM_RcppRMSE(SEXP vec1SEXP, SEXP vec2SEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppRMSE(vec1, vec2, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppAbsDiff
Rcpp::NumericVector RcppAbsDiff(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2);
RcppExport SEXP _spEDM_RcppAbsDiff(SEXP vec1SEXP, SEXP vec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec2(vec2SEXP);
    rcpp_result_gen = Rcpp::wrap(RcppAbsDiff(vec1, vec2));
    return rcpp_result_gen;
END_RCPP
}
// RcppSumNormalize
Rcpp::NumericVector RcppSumNormalize(const Rcpp::NumericVector& vec, bool NA_rm);
RcppExport SEXP _spEDM_RcppSumNormalize(SEXP vecSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSumNormalize(vec, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppDistance
double RcppDistance(const Rcpp::NumericVector& vec1, const Rcpp::NumericVector& vec2, bool L1norm, bool NA_rm);
RcppExport SEXP _spEDM_RcppDistance(SEXP vec1SEXP, SEXP vec2SEXP, SEXP L1normSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< bool >::type L1norm(L1normSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppDistance(vec1, vec2, L1norm, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppPearsonCor
double RcppPearsonCor(const Rcpp::NumericVector& y, const Rcpp::NumericVector& y_hat, bool NA_rm);
RcppExport SEXP _spEDM_RcppPearsonCor(SEXP ySEXP, SEXP y_hatSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_hat(y_hatSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppPearsonCor(y, y_hat, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppPartialCor
double RcppPartialCor(const Rcpp::NumericVector& y, const Rcpp::NumericVector& y_hat, const Rcpp::NumericMatrix& controls, bool NA_rm, bool linear);
RcppExport SEXP _spEDM_RcppPartialCor(SEXP ySEXP, SEXP y_hatSEXP, SEXP controlsSEXP, SEXP NA_rmSEXP, SEXP linearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_hat(y_hatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type controls(controlsSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    Rcpp::traits::input_parameter< bool >::type linear(linearSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppPartialCor(y, y_hat, controls, NA_rm, linear));
    return rcpp_result_gen;
END_RCPP
}
// RcppPartialCorTrivar
double RcppPartialCorTrivar(const Rcpp::NumericVector& y, const Rcpp::NumericVector& y_hat, const Rcpp::NumericVector& control, bool NA_rm, bool linear);
RcppExport SEXP _spEDM_RcppPartialCorTrivar(SEXP ySEXP, SEXP y_hatSEXP, SEXP controlSEXP, SEXP NA_rmSEXP, SEXP linearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y_hat(y_hatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type control(controlSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    Rcpp::traits::input_parameter< bool >::type linear(linearSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppPartialCorTrivar(y, y_hat, control, NA_rm, linear));
    return rcpp_result_gen;
END_RCPP
}
// RcppCorSignificance
double RcppCorSignificance(double r, int n, int k);
RcppExport SEXP _spEDM_RcppCorSignificance(SEXP rSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppCorSignificance(r, n, k));
    return rcpp_result_gen;
END_RCPP
}
// RcppCorConfidence
Rcpp::NumericVector RcppCorConfidence(double r, int n, int k, double level);
RcppExport SEXP _spEDM_RcppCorConfidence(SEXP rSEXP, SEXP nSEXP, SEXP kSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppCorConfidence(r, n, k, level));
    return rcpp_result_gen;
END_RCPP
}
// RcppMatDistance
Rcpp::NumericMatrix RcppMatDistance(const Rcpp::NumericMatrix& mat, bool L1norm, bool NA_rm);
RcppExport SEXP _spEDM_RcppMatDistance(SEXP matSEXP, SEXP L1normSEXP, SEXP NA_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< bool >::type L1norm(L1normSEXP);
    Rcpp::traits::input_parameter< bool >::type NA_rm(NA_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppMatDistance(mat, L1norm, NA_rm));
    return rcpp_result_gen;
END_RCPP
}
// RcppKNNIndice
Rcpp::IntegerVector RcppKNNIndice(const Rcpp::NumericMatrix& embedding_space, int target_idx, int k);
RcppExport SEXP _spEDM_RcppKNNIndice(SEXP embedding_spaceSEXP, SEXP target_idxSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type embedding_space(embedding_spaceSEXP);
    Rcpp::traits::input_parameter< int >::type target_idx(target_idxSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppKNNIndice(embedding_space, target_idx, k));
    return rcpp_result_gen;
END_RCPP
}
// RcppDistKNNIndice
Rcpp::IntegerVector RcppDistKNNIndice(const Rcpp::NumericMatrix& dist_mat, int target_idx, int k);
RcppExport SEXP _spEDM_RcppDistKNNIndice(SEXP dist_matSEXP, SEXP target_idxSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type dist_mat(dist_matSEXP);
    Rcpp::traits::input_parameter< int >::type target_idx(target_idxSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppDistKNNIndice(dist_mat, target_idx, k));
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
// RcppSVD
Rcpp::List RcppSVD(const Rcpp::NumericMatrix& X);
RcppExport SEXP _spEDM_RcppSVD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppSVD(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spEDM_RcppSimplexForecast", (DL_FUNC) &_spEDM_RcppSimplexForecast, 5},
    {"_spEDM_RcppSMapForecast", (DL_FUNC) &_spEDM_RcppSMapForecast, 6},
    {"_spEDM_RcppLaggedVar4Grid", (DL_FUNC) &_spEDM_RcppLaggedVar4Grid, 2},
    {"_spEDM_RcppGenGridEmbeddings", (DL_FUNC) &_spEDM_RcppGenGridEmbeddings, 3},
    {"_spEDM_RcppLocateGridIndices", (DL_FUNC) &_spEDM_RcppLocateGridIndices, 4},
    {"_spEDM_RcppSimplex4Grid", (DL_FUNC) &_spEDM_RcppSimplex4Grid, 7},
    {"_spEDM_RcppSMap4Grid", (DL_FUNC) &_spEDM_RcppSMap4Grid, 8},
    {"_spEDM_RcppGCCM4Grid", (DL_FUNC) &_spEDM_RcppGCCM4Grid, 11},
    {"_spEDM_RcppSCPCM4Grid", (DL_FUNC) &_spEDM_RcppSCPCM4Grid, 13},
    {"_spEDM_RcppGCMC4Grid", (DL_FUNC) &_spEDM_RcppGCMC4Grid, 9},
    {"_spEDM_DetectMaxNumThreads", (DL_FUNC) &_spEDM_DetectMaxNumThreads, 0},
    {"_spEDM_OptEmdedDim", (DL_FUNC) &_spEDM_OptEmdedDim, 1},
    {"_spEDM_OptThetaParm", (DL_FUNC) &_spEDM_OptThetaParm, 1},
    {"_spEDM_RcppLaggedNeighbor4Lattice", (DL_FUNC) &_spEDM_RcppLaggedNeighbor4Lattice, 2},
    {"_spEDM_RcppLaggedVar4Lattice", (DL_FUNC) &_spEDM_RcppLaggedVar4Lattice, 2},
    {"_spEDM_RcppGenLatticeEmbeddings", (DL_FUNC) &_spEDM_RcppGenLatticeEmbeddings, 4},
    {"_spEDM_RcppSimplex4Lattice", (DL_FUNC) &_spEDM_RcppSimplex4Lattice, 8},
    {"_spEDM_RcppSMap4Lattice", (DL_FUNC) &_spEDM_RcppSMap4Lattice, 9},
    {"_spEDM_RcppGCCM4Lattice", (DL_FUNC) &_spEDM_RcppGCCM4Lattice, 13},
    {"_spEDM_RcppSCPCM4Lattice", (DL_FUNC) &_spEDM_RcppSCPCM4Lattice, 15},
    {"_spEDM_RcppGCMC4Lattice", (DL_FUNC) &_spEDM_RcppGCMC4Lattice, 10},
    {"_spEDM_RcppMean", (DL_FUNC) &_spEDM_RcppMean, 2},
    {"_spEDM_RcppSum", (DL_FUNC) &_spEDM_RcppSum, 2},
    {"_spEDM_RcppVariance", (DL_FUNC) &_spEDM_RcppVariance, 2},
    {"_spEDM_RcppCovariance", (DL_FUNC) &_spEDM_RcppCovariance, 3},
    {"_spEDM_RcppMAE", (DL_FUNC) &_spEDM_RcppMAE, 3},
    {"_spEDM_RcppRMSE", (DL_FUNC) &_spEDM_RcppRMSE, 3},
    {"_spEDM_RcppAbsDiff", (DL_FUNC) &_spEDM_RcppAbsDiff, 2},
    {"_spEDM_RcppSumNormalize", (DL_FUNC) &_spEDM_RcppSumNormalize, 2},
    {"_spEDM_RcppDistance", (DL_FUNC) &_spEDM_RcppDistance, 4},
    {"_spEDM_RcppPearsonCor", (DL_FUNC) &_spEDM_RcppPearsonCor, 3},
    {"_spEDM_RcppPartialCor", (DL_FUNC) &_spEDM_RcppPartialCor, 5},
    {"_spEDM_RcppPartialCorTrivar", (DL_FUNC) &_spEDM_RcppPartialCorTrivar, 5},
    {"_spEDM_RcppCorSignificance", (DL_FUNC) &_spEDM_RcppCorSignificance, 3},
    {"_spEDM_RcppCorConfidence", (DL_FUNC) &_spEDM_RcppCorConfidence, 4},
    {"_spEDM_RcppMatDistance", (DL_FUNC) &_spEDM_RcppMatDistance, 3},
    {"_spEDM_RcppKNNIndice", (DL_FUNC) &_spEDM_RcppKNNIndice, 3},
    {"_spEDM_RcppDistKNNIndice", (DL_FUNC) &_spEDM_RcppDistKNNIndice, 3},
    {"_spEDM_RcppLinearTrendRM", (DL_FUNC) &_spEDM_RcppLinearTrendRM, 4},
    {"_spEDM_RcppSVD", (DL_FUNC) &_spEDM_RcppSVD, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_spEDM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

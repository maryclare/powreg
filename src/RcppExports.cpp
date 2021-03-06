// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rtnormrej
NumericVector rtnormrej(NumericVector mu, NumericVector sd, NumericVector l, NumericVector r);
RcppExport SEXP _powreg_rtnormrej(SEXP muSEXP, SEXP sdSEXP, SEXP lSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(rtnormrej(mu, sd, l, r));
    return rcpp_result_gen;
END_RCPP
}
// rshiftexp
NumericVector rshiftexp(NumericVector d, NumericVector t);
RcppExport SEXP _powreg_rshiftexp(SEXP dSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(rshiftexp(d, t));
    return rcpp_result_gen;
END_RCPP
}
// remcol
arma::mat remcol(arma::mat A, int i);
RcppExport SEXP _powreg_remcol(SEXP ASEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(remcol(A, i));
    return rcpp_result_gen;
END_RCPP
}
// remrow
arma::colvec remrow(arma::colvec a, int i);
RcppExport SEXP _powreg_remrow(SEXP aSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(remrow(a, i));
    return rcpp_result_gen;
END_RCPP
}
// sampleBeta
arma::colvec sampleBeta(NumericVector start, NumericVector DUty, NumericVector delta, NumericVector d, NumericMatrix Vt, double sigsq, NumericMatrix W);
RcppExport SEXP _powreg_sampleBeta(SEXP startSEXP, SEXP DUtySEXP, SEXP deltaSEXP, SEXP dSEXP, SEXP VtSEXP, SEXP sigsqSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type DUty(DUtySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Vt(VtSEXP);
    Rcpp::traits::input_parameter< double >::type sigsq(sigsqSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleBeta(start, DUty, delta, d, Vt, sigsq, W));
    return rcpp_result_gen;
END_RCPP
}
// sampleGamma
NumericVector sampleGamma(NumericVector beta, double tausq, double q);
RcppExport SEXP _powreg_sampleGamma(SEXP betaSEXP, SEXP tausqSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type tausq(tausqSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleGamma(beta, tausq, q));
    return rcpp_result_gen;
END_RCPP
}
// sampler
List sampler(const NumericVector& DUty, const NumericMatrix& Vt, const NumericVector& d, const NumericMatrix& W, double sigsq, double tausq, double q, const int& samples, NumericVector start, int seed, const int& burn, const int& thin);
RcppExport SEXP _powreg_sampler(SEXP DUtySEXP, SEXP VtSEXP, SEXP dSEXP, SEXP WSEXP, SEXP sigsqSEXP, SEXP tausqSEXP, SEXP qSEXP, SEXP samplesSEXP, SEXP startSEXP, SEXP seedSEXP, SEXP burnSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type DUty(DUtySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Vt(VtSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type sigsq(sigsqSEXP);
    Rcpp::traits::input_parameter< double >::type tausq(tausqSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler(DUty, Vt, d, W, sigsq, tausq, q, samples, start, seed, burn, thin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_powreg_rtnormrej", (DL_FUNC) &_powreg_rtnormrej, 4},
    {"_powreg_rshiftexp", (DL_FUNC) &_powreg_rshiftexp, 2},
    {"_powreg_remcol", (DL_FUNC) &_powreg_remcol, 2},
    {"_powreg_remrow", (DL_FUNC) &_powreg_remrow, 2},
    {"_powreg_sampleBeta", (DL_FUNC) &_powreg_sampleBeta, 7},
    {"_powreg_sampleGamma", (DL_FUNC) &_powreg_sampleGamma, 3},
    {"_powreg_sampler", (DL_FUNC) &_powreg_sampler, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_powreg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

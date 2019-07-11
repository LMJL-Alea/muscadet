// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEnsmallen.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EvaluateBessel
double EvaluateBessel(const arma::vec& p, const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double amplitude1, const double amplitude2, const double amplitude12, const double alpha1, const double alpha2, const double alpha12);
RcppExport SEXP _mediator_EvaluateBessel(SEXP pSEXP, SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP amplitude1SEXP, SEXP amplitude2SEXP, SEXP amplitude12SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP alpha12SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude1(amplitude1SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude2(amplitude2SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude12(amplitude12SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha12(alpha12SEXP);
    rcpp_result_gen = Rcpp::wrap(EvaluateBessel(p, X, labels, lb, ub, amplitude1, amplitude2, amplitude12, alpha1, alpha2, alpha12));
    return rcpp_result_gen;
END_RCPP
}
// Estimate
arma::mat Estimate(const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double amplitude1, const double amplitude2, const double amplitude12, const double alpha1, const double alpha2, const double alpha12);
RcppExport SEXP _mediator_Estimate(SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP amplitude1SEXP, SEXP amplitude2SEXP, SEXP amplitude12SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP alpha12SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude1(amplitude1SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude2(amplitude2SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude12(amplitude12SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha12(alpha12SEXP);
    rcpp_result_gen = Rcpp::wrap(Estimate(X, labels, lb, ub, amplitude1, amplitude2, amplitude12, alpha1, alpha2, alpha12));
    return rcpp_result_gen;
END_RCPP
}
// Evaluate
double Evaluate(const arma::vec& p, const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double amplitude1, const double amplitude2, const double amplitude12, const double alpha1, const double alpha2, const double alpha12);
RcppExport SEXP _mediator_Evaluate(SEXP pSEXP, SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP amplitude1SEXP, SEXP amplitude2SEXP, SEXP amplitude12SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP alpha12SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude1(amplitude1SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude2(amplitude2SEXP);
    Rcpp::traits::input_parameter< const double >::type amplitude12(amplitude12SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha12(alpha12SEXP);
    rcpp_result_gen = Rcpp::wrap(Evaluate(p, X, labels, lb, ub, amplitude1, amplitude2, amplitude12, alpha1, alpha2, alpha12));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mediator_EvaluateBessel", (DL_FUNC) &_mediator_EvaluateBessel, 11},
    {"_mediator_Estimate", (DL_FUNC) &_mediator_Estimate, 10},
    {"_mediator_Evaluate", (DL_FUNC) &_mediator_Evaluate, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_mediator(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

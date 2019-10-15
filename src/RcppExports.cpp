// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEnsmallen.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EstimateBessel
arma::mat EstimateBessel(const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double rho2, const double alpha1, const double alpha2, const bool estimate_alpha);
RcppExport SEXP _mediator_EstimateBessel(SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP estimate_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type estimate_alpha(estimate_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(EstimateBessel(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha));
    return rcpp_result_gen;
END_RCPP
}
// EvaluateBessel
double EvaluateBessel(const arma::vec& p, const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double alpha1, const double rho2, const double alpha2);
RcppExport SEXP _mediator_EvaluateBessel(SEXP pSEXP, SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP alpha1SEXP, SEXP rho2SEXP, SEXP alpha2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    rcpp_result_gen = Rcpp::wrap(EvaluateBessel(p, X, labels, lb, ub, rho1, alpha1, rho2, alpha2));
    return rcpp_result_gen;
END_RCPP
}
// InitializeBessel
arma::mat InitializeBessel(const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double rho2, const double alpha1, const double alpha2, const bool estimate_alpha);
RcppExport SEXP _mediator_InitializeBessel(SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP estimate_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type estimate_alpha(estimate_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(InitializeBessel(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha));
    return rcpp_result_gen;
END_RCPP
}
// EstimateGauss
arma::mat EstimateGauss(const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double rho2, const double alpha1, const double alpha2, const bool estimate_alpha);
RcppExport SEXP _mediator_EstimateGauss(SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP estimate_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type estimate_alpha(estimate_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(EstimateGauss(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha));
    return rcpp_result_gen;
END_RCPP
}
// EvaluateGauss
double EvaluateGauss(const arma::vec& p, const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double rho2, const double alpha1, const double alpha2, const bool estimate_alpha);
RcppExport SEXP _mediator_EvaluateGauss(SEXP pSEXP, SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP estimate_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type estimate_alpha(estimate_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(EvaluateGauss(p, X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha));
    return rcpp_result_gen;
END_RCPP
}
// InitializeGauss
arma::mat InitializeGauss(const arma::mat& X, const arma::uvec& labels, const arma::vec& lb, const arma::vec& ub, const double rho1, const double rho2, const double alpha1, const double alpha2, const bool estimate_alpha);
RcppExport SEXP _mediator_InitializeGauss(SEXP XSEXP, SEXP labelsSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP estimate_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< const double >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< const double >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type estimate_alpha(estimate_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(InitializeGauss(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mediator_EstimateBessel", (DL_FUNC) &_mediator_EstimateBessel, 9},
    {"_mediator_EvaluateBessel", (DL_FUNC) &_mediator_EvaluateBessel, 9},
    {"_mediator_InitializeBessel", (DL_FUNC) &_mediator_InitializeBessel, 9},
    {"_mediator_EstimateGauss", (DL_FUNC) &_mediator_EstimateGauss, 9},
    {"_mediator_EvaluateGauss", (DL_FUNC) &_mediator_EvaluateGauss, 10},
    {"_mediator_InitializeGauss", (DL_FUNC) &_mediator_InitializeGauss, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_mediator(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

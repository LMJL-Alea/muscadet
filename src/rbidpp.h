#include <RcppArmadillo.h>

arma::mat22 GetGaussianKernel(
    const double rSq,
    const double rho1,
    const double rho2,
    const double alpha1Sq,
    const double alpha2Sq,
    const double alpha12Sq,
    const double tau
);

arma::mat22 GetBesselKernel(
    const double rSq,
    const double rho1,
    const double rho2,
    const double alpha1Sq,
    const double alpha2Sq,
    const double alpha12Sq,
    const double tau
);

// [[Rcpp::export]]
Rcpp::List rbidpp_impl(
    const int N,
    const double L,
    const double rho1,
    const double rho2,
    const double alpha1,
    const double alpha2,
    const double alpha12,
    const double tau,
    const std::string model,
    const unsigned int nbThreads
);

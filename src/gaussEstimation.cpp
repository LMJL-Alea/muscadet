#include <RcppEnsmallen.h>
#include "logLikelihood.h"

//' Stationary Bivariate Gaussian DPP Estimator
//'
//' This function estimates the parameters of a stationary bivariate Gaussian DPP from a set of observed points and labels.
//'
//' @param X A matrix of size n x (d+1) storing the points in R^d and their label in last column.
//'
//' @return A vector with the estimated model parameters.
//'
//' @export
//' @examples
//' # Simulate some data
//' m <- spatstat::dppGauss(lambda = 100, alpha = 0.05, d = 2)
//' X1 <- stats::simulate(m)
//' X2 <- stats::simulate(m)
//' X1 <- cbind(X1$x, X1$y, rep(1, X1$n))
//' X2 <- cbind(X2$x, X2$y, rep(2, X2$n))
//' X <- rbind(X1, X2)
//'
//' # Run parameter estimation
//' params <- Estimate(X)
//'
//' # Verify parameters were recovered
//' params
// [[Rcpp::export]]
arma::mat Estimate(const arma::mat& X) {

  // Construct the objective function.
  GaussianLogLikelihood logLik;
  logLik.SetInputs(X);

  // Create the Augmented Lagrangian optimizer with default parameters.
  // The ens::L_BFGS is used internally.
  ens::AugLagrangian auglag;

  // Create a starting point for our optimization randomly.
  // The model has p parameters, so the shape is p x 1.
  // arma::mat params(4, 1, arma::fill::randn);
  arma::mat params(4, 1);
  params[0] = std::log(0.05);
  params[1] = std::log(0.05);
  params[2] = std::log(0.05);
  params[3] = 0.001;

  // Time the routine
  arma::wall_clock clock;
  clock.tic();

  // Run the optimization
  Rcpp::Rcout << "Initial parameters: " << params << std::endl;
  auglag.Optimize(logLik, params);
  Rcpp::Rcout << "Final parameters: " << params << std::endl;

  // End time output
  Rcpp::Rcout << "Estimation performed in " << clock.toc() << " seconds." << std::endl;

  for (unsigned int i = 0;i < 3;++i)
    params[i] = std::exp(params[i]);

  return params;
}

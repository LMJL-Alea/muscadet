#include <RcppEnsmallen.h>
#include "gaussianLogLikelihood.h"

arma::mat GetStartingPoint(const double rho1, const double rho2, const double alpha1, const double alpha2, const unsigned int dimension)
{
  const double epsilon = 1.0e-4;
  arma::mat params(4, 1);

  double alpha12 = 2.0 * std::sqrt((alpha1 * alpha1 + alpha2 * alpha2) / 2.0);

  double amp1 = rho1 * std::pow(std::sqrt(M_PI) * alpha1, (double)dimension);
  double amp2 = rho2 * std::pow(std::sqrt(M_PI) * alpha2, (double)dimension);
  double amp12 = std::pow(std::sqrt(M_PI) * alpha12, (double)dimension);
  double ub = std::sqrt(amp1 * amp2 / amp12);
  ub = std::min(ub, std::sqrt(4.0 * (1.0 - amp1) * (1.0 - amp2) / amp12)) - epsilon;

  double lb = -ub;
  double sigma12 = lb + ub * arma::randu();

  params[0] = std::log(alpha1);
  params[1] = std::log(alpha12);
  params[2] = std::log(alpha2);
  params[3] = sigma12;

  return params;
}

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
//' library(spatstat)
//' # Simulate some data
//' m <- dppGauss(lambda = 100, alpha = 0.05, d = 2)
//' X1 <- simulate(m)
//' fit1 <- dppm(X1, dppGauss, method = "palm")$fitted$fixedpar
//' X2 <- stats::simulate(m)
//' fit2 <- dppm(X2, dppGauss, method = "palm")$fitted$fixedpar
//' X1 <- cbind(X1$x, X1$y, rep(1, X1$n))
//' X2 <- cbind(X2$x, X2$y, rep(2, X2$n))
//' X <- rbind(X1, X2)
//'
//' # Run parameter estimation
//' params <- Estimate(X, fit1$lambda, fit2$lambda, fit1$alpha, fit2$alpha)
//'
//' # Verify parameters were recovered
//' params
// [[Rcpp::export]]
arma::mat Estimate(const arma::mat& X, const double rho1, const double rho2, const double alpha1, const double alpha2, const double volume = 1.0)
{
  // Construct the objective function.
  GaussianLogLikelihood logLik;
  logLik.SetInputs(X, rho1, rho2, volume);

  // Create the Augmented Lagrangian optimizer with default parameters.
  // The ens::L_BFGS is used internally.
  // ens::AugLagrangian optimizer;
  // ens::DE optimizer;
  // ens::CNE optimizer;
  ens::L_BFGS optimizer;

  // Create a starting point for our optimization randomly within the
  // authorized search space.

  arma::mat params = GetStartingPoint(rho1, rho2, alpha1, alpha2, X.n_cols - 1);

  // Time the routine
  arma::wall_clock clock;
  clock.tic();

  // Run the optimization
  Rcpp::Rcout << "Initial parameters: " << params.as_row() << std::endl;
  optimizer.Optimize(logLik, params);
  Rcpp::Rcout << "Final parameters: " << params.as_row() << std::endl;

  // End time output
  Rcpp::Rcout << "Estimation performed in " << clock.toc() << " seconds." << std::endl;

  for (unsigned int i = 0;i < 3;++i)
    params[i] = std::exp(params[i]);

  return params;
}

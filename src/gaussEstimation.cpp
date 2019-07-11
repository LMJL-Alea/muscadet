#include <RcppEnsmallen.h>
#include "gaussianLogLikelihood.h"

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
arma::mat Estimate(
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1,
    const double rho2,
    const double alpha1,
    const double alpha2)
{
  // Construct the objective function.
  // GaussianLogLikelihood logLik;
  GaussianLogLikelihoodV2 logLik;
  logLik.SetInputs(X, labels, lb, ub);
  logLik.SetFirstIntensity(rho1);
  logLik.SetSecondIntensity(rho2);
  logLik.SetFirstAlpha(alpha1);
  logLik.SetSecondAlpha(alpha2);
  logLik.SetCrossAlpha(0.03);

  // Create the Augmented Lagrangian optimizer with default parameters.
  // The ens::L_BFGS is used internally.
  // ens::AugLagrangian optimizer;
  ens::DE optimizer;
  // ens::SPSA optimizer(0.1, 0.102, 0.16, 0.3, 100000, 1e-5);
  // ens::ExponentialSchedule expSchedule;
  // ens::SA<> optimizer(expSchedule);
  // ens::CNE optimizer(200, 10000, 0.2, 0.2, 0.3, 1e-5);
  // ens::L_BFGS optimizer;

  // Create a starting point for our optimization randomly within the
  // authorized search space.
  // arma::mat params = logLik.GetInitialPoint(rho1, rho2, alpha1, alpha2);
  arma::mat params(1, 1);
  // params[0] = 1.0001 * std::sqrt((alpha1 * alpha1 + alpha2 * alpha2) / 2.0);
  // params[1] = 0.7;
  params[0] = 0.2;

  // Time the routine
  arma::wall_clock clock;
  clock.tic();

  // Run the optimization
  Rcpp::Rcout << "Initial parameters: " << params.as_row() << std::endl;
  optimizer.Optimize(logLik, params);
  Rcpp::Rcout << "Final parameters: " << params.as_row() << std::endl;

  // End time output
  Rcpp::Rcout << "Estimation performed in " << clock.toc() << " seconds." << std::endl;
  Rcpp::Rcout << "Min: " << logLik.Evaluate(params) << std::endl;

  // for (unsigned int i = 0;i < 3;++i)
  //   params[i] = std::exp(params[i]);

  return params;
}

// [[Rcpp::export]]
double Evaluate(
    const arma::vec &p,
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1,
    const double rho2,
    const double alpha1,
    const double alpha2)
{
  // Construct the objective function.
  GaussianLogLikelihoodV2 logLik;
  logLik.SetInputs(X, labels, lb, ub);

  logLik.SetFirstAlpha(alpha1);
  logLik.SetFirstIntensity(rho1);

  logLik.SetSecondAlpha(alpha2);
  logLik.SetSecondIntensity(rho2);

  // logLik.SetCrossAlpha(0.04);
  logLik.SetCrossIntensity(50.0 * M_PI * 0.04 * 0.04);

  arma::mat params(p.n_elem, 1);
  for (unsigned int i = 0;i < p.n_elem;++i)
    params[i] = p[i];
  return logLik.Evaluate(params);
}

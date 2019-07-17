#include <RcppEnsmallen.h>
#include "besselLogLikelihood.h"

//' Stationary Bivariate Bessel DPP Estimator
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
arma::mat EstimateBessel(
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double amplitude1 = NA_REAL,
    const double amplitude2 = NA_REAL,
    const double amplitude12 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const double alpha12 = NA_REAL)
{
  // Construct the objective function.
  BesselLogLikelihood logLik;

  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
    logLik.SetFirstAlpha(alpha1);

  if (arma::is_finite(amplitude1))
    logLik.SetFirstAmplitude(amplitude1);

  if (arma::is_finite(alpha2))
    logLik.SetSecondAlpha(alpha2);

  if (arma::is_finite(amplitude2))
    logLik.SetSecondAmplitude(amplitude2);

  if (arma::is_finite(alpha12))
    logLik.SetCrossAlpha(alpha12);

  if (arma::is_finite(amplitude12))
    logLik.SetCrossAmplitude(amplitude12);

  // Create the Augmented Lagrangian optimizer with default parameters.
  // The ens::L_BFGS is used internally.
  // ens::AugLagrangian optimizer;
  // ens::DE optimizer;
  // ens::SPSA optimizer(0.1, 0.102, 0.16, 0.3, 100000, 1e-5);
  ens::ExponentialSchedule expSchedule;
  ens::SA<> optimizer(expSchedule);
  // ens::CNE optimizer(200, 10000, 0.2, 0.2, 0.3, 1e-5);
  // ens::L_BFGS optimizer;

  // Create a starting point for our optimization randomly within the
  // authorized search space.
  arma::mat params = logLik.GetInitialPoint();

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

  return params;
}

// [[Rcpp::export]]
double EvaluateBessel(
    const arma::vec &p,
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double amplitude1 = NA_REAL,
    const double amplitude2 = NA_REAL,
    const double amplitude12 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const double alpha12 = NA_REAL)
{
  // Construct the objective function.
  BesselLogLikelihood logLik;

  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
    logLik.SetFirstAlpha(alpha1);

  if (arma::is_finite(amplitude1))
    logLik.SetFirstAmplitude(amplitude1);

  if (arma::is_finite(alpha2))
    logLik.SetSecondAlpha(alpha2);

  if (arma::is_finite(amplitude2))
    logLik.SetSecondAmplitude(amplitude2);

  if (arma::is_finite(alpha12))
    logLik.SetCrossAlpha(alpha12);

  if (arma::is_finite(amplitude12))
    logLik.SetCrossAmplitude(amplitude12);

  arma::mat params(p.n_elem, 1);
  for (unsigned int i = 0;i < p.n_elem;++i)
    params[i] = p[i];

  return logLik.Evaluate(params);
}

// [[Rcpp::export]]
arma::mat InitializeBessel(
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double amplitude1 = NA_REAL,
    const double amplitude2 = NA_REAL,
    const double amplitude12 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const double alpha12 = NA_REAL)
{
  // Construct the objective function.
  BesselLogLikelihood logLik;
  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
    logLik.SetFirstAlpha(alpha1);

  if (arma::is_finite(amplitude1))
    logLik.SetFirstAmplitude(amplitude1);

  if (arma::is_finite(alpha2))
    logLik.SetSecondAlpha(alpha2);

  if (arma::is_finite(amplitude2))
    logLik.SetSecondAmplitude(amplitude2);

  if (arma::is_finite(alpha12))
    logLik.SetCrossAlpha(alpha12);

  if (arma::is_finite(amplitude12))
    logLik.SetCrossAmplitude(amplitude12);

  return logLik.GetInitialPoint();
}

#include <RcppEnsmallen.h>
#include "gaussLogLikelihood.h"

//' Stationary Bivariate Gaussian DPP Estimator
//'
//' This function estimates the parameters of a stationary bivariate Gaussian DPP from a set of observed points and labels.
//'
//' @param X A matrix of size n x (d+1) storing the points in R^d and their label in last column.
//'
//' @return A vector with the estimated model parameters.
//'
//' @keywords internal
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
arma::mat EstimateGauss(
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1 = NA_REAL,
    const double rho2 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const bool estimate_alpha = true)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;
  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
  {
    if (!estimate_alpha)
      logLik.SetFirstAlpha(alpha1);

    if (arma::is_finite(rho1))
    {
      double a1 = logLik.RetrieveAmplitudeFromParameters(rho1, alpha1, X.n_cols);
      // logLik.SetFirstAmplitude(a1);
    }
  }

  if (arma::is_finite(alpha2))
  {
    if (!estimate_alpha)
      logLik.SetSecondAlpha(alpha2);

    if (arma::is_finite(rho2))
    {
      double a2 = logLik.RetrieveAmplitudeFromParameters(rho2, alpha2, X.n_cols);
      // logLik.SetSecondAmplitude(a2);
    }
  }

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
double EvaluateGauss(
    const arma::vec &p,
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1 = NA_REAL,
    const double rho2 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const bool estimate_alpha = true)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;
  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
  {
    if (!estimate_alpha)
      logLik.SetFirstAlpha(alpha1);

    if (arma::is_finite(rho1))
    {
      double a1 = logLik.RetrieveAmplitudeFromParameters(rho1, alpha1, X.n_cols);
      // logLik.SetFirstAmplitude(a1);
    }
  }

  if (arma::is_finite(alpha2))
  {
    if (!estimate_alpha)
      logLik.SetSecondAlpha(alpha2);

    if (arma::is_finite(rho2))
    {
      double a2 = logLik.RetrieveAmplitudeFromParameters(rho2, alpha2, X.n_cols);
      // logLik.SetSecondAmplitude(a2);
    }
  }

  arma::mat params(p.n_elem, 1);
  for (unsigned int i = 0;i < p.n_elem;++i)
    params[i] = p[i];

  return logLik.Evaluate(params);
}

// [[Rcpp::export]]
arma::mat InitializeGauss(
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1 = NA_REAL,
    const double rho2 = NA_REAL,
    const double alpha1 = NA_REAL,
    const double alpha2 = NA_REAL,
    const bool estimate_alpha = true)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;
  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(alpha1))
  {
    if (!estimate_alpha)
      logLik.SetFirstAlpha(alpha1);

    if (arma::is_finite(rho1))
    {
      double a1 = logLik.RetrieveAmplitudeFromParameters(rho1, alpha1, X.n_cols);
      // logLik.SetFirstAmplitude(a1);
    }
  }

  if (arma::is_finite(alpha2))
  {
    if (!estimate_alpha)
      logLik.SetSecondAlpha(alpha2);

    if (arma::is_finite(rho2))
    {
      double a2 = logLik.RetrieveAmplitudeFromParameters(rho2, alpha2, X.n_cols);
      // logLik.SetSecondAmplitude(a2);
    }
  }

  return logLik.GetInitialPoint();
}

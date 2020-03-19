#include <RcppEnsmallen.h>
#include "gaussLogLikelihood.h"

//' Stationary Bivariate Gaussian DPP Estimator
//'
//' This function estimates the parameters of a stationary bivariate Gaussian
//' DPP from a set of observed points and labels.
//'
//' @param X A matrix of size n x (d+1) storing the points in R^d and their label in last column.
//'
//' @return A vector with the estimated model parameters.
//'
//' @export
//' @examples
//' dpp <- sim_gauss5[[1]]
//' X <- cbind(dpp$x, dpp$y)
//' labels <- dpp$marks
//' rho1 <- rho2 <- 100
//' rho12 <- sqrt(0.5 * sqrt(rho1 * rho2)
//' alpha1 <- alpha2 <- 0.03
//' alpha12 <- 0.03
//' d <- 2
//' EstimateGauss(
//'   X = X,
//'   labels = labels,
//'   lb = rep(-0.5, ncol(X)),
//'   ub = rep( 0.5, ncol(X)),
//'   rho1 = rho1,
//'   rho2 = rho2,
//'   alpha1 = alpha1,
//'   alpha2 = alpha2,
//'   estimate_alpha = FALSE
//' )
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

  // Create the optimizer.
  // ens::DE optimizer(1000, 1000, 0.6, 0.8, 1e-5);
  ens::ExponentialSchedule expSchedule;
  ens::SA<> optimizer(expSchedule);//, 1000000, 1000., 1000, 100, 1e-10, 3, 1.5, 0.5, 0.3);
  // ens::CNE optimizer(200, 10000, 0.2, 0.2, 0.3, 1e-5);

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

//' @export
// [[Rcpp::export]]
double EvaluateGauss(
    const arma::vec &p,
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1 = NA_REAL,
    const double rho2 = NA_REAL)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;
  logLik.SetInputs(X, labels, lb, ub);

  if (arma::is_finite(rho1) & arma::is_finite(rho2))
    logLik.SetIntensities(rho1, rho2);

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
  return logLik.GetInitialPoint();
}

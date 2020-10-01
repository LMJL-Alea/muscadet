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
    const arma::mat &points,
    const double alpha1,
    const Rcpp::Nullable<double> &rho1 = R_NilValue,
    const Rcpp::Nullable<double> &alpha2 = R_NilValue,
    const Rcpp::Nullable<double> &rho2 = R_NilValue,
    const Rcpp::Nullable<double> &alpha12 = R_NilValue,
    const Rcpp::Nullable<double> &tau = R_NilValue,
    const Rcpp::Nullable<arma::vec> &lower_bound = R_NilValue,
    const Rcpp::Nullable<arma::vec> &upper_bound = R_NilValue,
    const Rcpp::Nullable<arma::uvec> &labels = R_NilValue,
    const int N = 50,
    const bool estimate_intensities = false,
    const bool use_verbose = false)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;

  arma::vec lb, ub;

  if (lower_bound.isNull())
  {
    lb.set_size(points.n_cols);
    lb.fill(0.0);
  }
  else
    lb = Rcpp::as<arma::vec>(lower_bound);

  if (upper_bound.isNull())
  {
    ub.set_size(points.n_cols);
    ub.fill(1.0);
  }
  else
    ub = Rcpp::as<arma::vec>(upper_bound);

  arma::uvec labs;

  if (labels.isNull())
  {
    labs.set_size(points.n_rows);
    labs.fill(0.0);
  }

  logLik.SetInputData(points, lb, ub, labs);
  logLik.SetTruncationIndex(N);
  logLik.SetUseVerbose(use_verbose);

  logLik.SetFirstAlpha(alpha1);

  if (alpha2.isNotNull())
    logLik.SetSecondAlpha(Rcpp::as<double>(alpha2));

  if (alpha12.isNotNull())
    logLik.SetCrossAlpha(Rcpp::as<double>(alpha12));

  if (tau.isNotNull())
    logLik.SetCorrelation(Rcpp::as<double>(tau));

  arma::mat params = logLik.GetInitialPoint();

  // Create the optimizer.
  ens::SPSA optimizer;
  // ens::ExponentialSchedule expSchedule;
  // ens::SA<> optimizer(expSchedule);//, 1000000, 1000., 1000, 100, 1e-10, 3, 1.5, 0.5, 0.3);
  // ens::CNE optimizer(200, 10000, 0.2, 0.2, 0.3, 1e-5);

  // Time the routine
  arma::wall_clock clock;
  clock.tic();

  // Run the optimization
  Rcpp::Rcout << "* Initial parameters: " << params.as_row() << std::endl;
  optimizer.Optimize(logLik, params);
  Rcpp::Rcout << "* Final parameters: " << params.as_row() << std::endl;

  // End time output
  Rcpp::Rcout << "Estimation performed in " << clock.toc() << " seconds." << std::endl;
  Rcpp::Rcout << "Min: " << logLik.Evaluate(params) << std::endl;

  return params;
}

//' @export
// [[Rcpp::export]]
double EvaluateGauss(
    const arma::mat &points,
    const double alpha1,
    const Rcpp::Nullable<double> &rho1 = R_NilValue,
    const Rcpp::Nullable<double> &alpha2 = R_NilValue,
    const Rcpp::Nullable<double> &rho2 = R_NilValue,
    const Rcpp::Nullable<double> &alpha12 = R_NilValue,
    const Rcpp::Nullable<double> &tau = R_NilValue,
    const Rcpp::Nullable<arma::vec> &lower_bound = R_NilValue,
    const Rcpp::Nullable<arma::vec> &upper_bound = R_NilValue,
    const Rcpp::Nullable<arma::uvec> &labels = R_NilValue,
    const int N = 50,
    const bool estimate_intensities = false)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;

  arma::vec lb, ub;

  if (lower_bound.isNull())
  {
    lb.set_size(points.n_cols);
    lb.fill(0.0);
  }
  else
    lb = Rcpp::as<arma::vec>(lower_bound);

  if (upper_bound.isNull())
  {
    ub.set_size(points.n_cols);
    ub.fill(1.0);
  }
  else
    ub = Rcpp::as<arma::vec>(upper_bound);

  arma::uvec labs;

  if (labels.isNull())
  {
    labs.set_size(points.n_rows);
    labs.fill(0.0);
  }

  logLik.SetInputData(points, lb, ub, labs);
  logLik.SetTruncationIndex(N);

  logLik.SetFirstAlpha(alpha1);

  if (alpha2.isNotNull())
    logLik.SetSecondAlpha(Rcpp::as<double>(alpha2));

  if (alpha12.isNotNull())
    logLik.SetCrossAlpha(Rcpp::as<double>(alpha12));

  if (tau.isNotNull())
    logLik.SetCorrelation(Rcpp::as<double>(tau));

  arma::mat params = logLik.GetInitialPoint();

  return logLik.Evaluate(params);
}

//' @export
// [[Rcpp::export]]
arma::mat InitializeGauss(const arma::mat &points,
                          const double alpha1,
                          const Rcpp::Nullable<double> &rho1 = R_NilValue,
                          const Rcpp::Nullable<double> &alpha2 = R_NilValue,
                          const Rcpp::Nullable<double> &rho2 = R_NilValue,
                          const Rcpp::Nullable<double> &alpha12 = R_NilValue,
                          const Rcpp::Nullable<double> &tau = R_NilValue,
                          const Rcpp::Nullable<arma::vec> &lower_bound = R_NilValue,
                          const Rcpp::Nullable<arma::vec> &upper_bound = R_NilValue,
                          const Rcpp::Nullable<arma::uvec> &labels = R_NilValue,
                          const int N = 50,
                          const bool estimate_intensities = false,
                          const bool use_verbose = false)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;

  arma::vec lb, ub;

  if (lower_bound.isNull())
  {
    lb.set_size(points.n_cols);
    lb.fill(0.0);
  }
  else
    lb = Rcpp::as<arma::vec>(lower_bound);

  if (upper_bound.isNull())
  {
    ub.set_size(points.n_cols);
    ub.fill(1.0);
  }
  else
    ub = Rcpp::as<arma::vec>(upper_bound);

  arma::uvec labs;

  if (labels.isNull())
  {
    labs.set_size(points.n_rows);
    labs.fill(0.0);
  }

  logLik.SetInputData(points, lb, ub, labs);
  logLik.SetTruncationIndex(N);
  logLik.SetUseVerbose(use_verbose);

  logLik.SetFirstAlpha(alpha1);

  if (alpha2.isNotNull())
    logLik.SetSecondAlpha(Rcpp::as<double>(alpha2));

  if (alpha12.isNotNull())
    logLik.SetCrossAlpha(Rcpp::as<double>(alpha12));

  if (tau.isNotNull())
    logLik.SetCorrelation(Rcpp::as<double>(tau));

  return logLik.GetInitialPoint();
}

//' @export
// [[Rcpp::export]]
double MLEGauss(const arma::vec &x,
                const arma::mat &points,
                const double alpha1,
                const Rcpp::Nullable<double> &rho1 = R_NilValue,
                const Rcpp::Nullable<double> &alpha2 = R_NilValue,
                const Rcpp::Nullable<double> &rho2 = R_NilValue,
                const Rcpp::Nullable<double> &alpha12 = R_NilValue,
                const Rcpp::Nullable<double> &tau = R_NilValue,
                const Rcpp::Nullable<arma::vec> &lower_bound = R_NilValue,
                const Rcpp::Nullable<arma::vec> &upper_bound = R_NilValue,
                const Rcpp::Nullable<arma::uvec> &labels = R_NilValue,
                const int N = 50,
                const bool estimate_intensities = false,
                const bool use_verbose = false)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;

  arma::vec lb, ub;

  if (lower_bound.isNull())
  {
    lb.set_size(points.n_cols);
    lb.fill(0.0);
  }
  else
    lb = Rcpp::as<arma::vec>(lower_bound);

  if (upper_bound.isNull())
  {
    ub.set_size(points.n_cols);
    ub.fill(1.0);
  }
  else
    ub = Rcpp::as<arma::vec>(upper_bound);

  arma::uvec labs;

  if (labels.isNull())
  {
    labs.set_size(points.n_rows);
    labs.fill(0.0);
  }

  logLik.SetInputData(points, lb, ub, labs);
  logLik.SetTruncationIndex(N);
  logLik.SetUseVerbose(use_verbose);

  logLik.SetFirstAlpha(alpha1);

  if (alpha2.isNotNull())
    logLik.SetSecondAlpha(Rcpp::as<double>(alpha2));

  if (alpha12.isNotNull())
    logLik.SetCrossAlpha(Rcpp::as<double>(alpha12));

  if (tau.isNotNull())
    logLik.SetCorrelation(Rcpp::as<double>(tau));

  return logLik.Evaluate(x);
}

//' @export
// [[Rcpp::export]]
arma::vec log_likelihood(const arma::mat &theta,
                         const arma::mat &points,
                         const arma::vec &lower_bound,
                         const arma::vec &upper_bound,
                         const Rcpp::Nullable<arma::uvec> &marks = R_NilValue,
                         const int N = 50,
                         const bool use_verbose = false)
{
  // Construct the objective function.
  GaussLogLikelihood logLik;

  arma::uvec pointMarks;

  if (marks.isNull())
  {
    pointMarks.set_size(points.n_rows);
    pointMarks.fill(1);
  }
  else
    pointMarks = Rcpp::as<arma::uvec>(marks);

  logLik.SetInputData(points, lower_bound, upper_bound, pointMarks);
  logLik.SetTruncationIndex(N);
  logLik.SetUseVerbose(use_verbose);

  unsigned int numParams = theta.n_rows;

  arma::vec outputValues(theta.n_cols);
  arma::mat workingParams(theta.n_rows, 1);

  arma::wall_clock timer;
  timer.tic();

  for (unsigned int i = 0;i < theta.n_cols;++i)
  {
    if (numParams == 1)
      logLik.SetFirstAlpha(theta(0, i));

    if (numParams == 4)
    {
      logLik.SetFirstAlpha(theta(0, i));
      logLik.SetSecondAlpha(theta(1, i));
      logLik.SetCrossAlpha(theta(2, i));
      logLik.SetCorrelation(theta(3, i));
    }

    workingParams = logLik.GetInitialPoint();
    outputValues(i) = logLik.Evaluate(workingParams);
  }

  Rcpp::Rcout << "It took " << timer.toc() << " seconds for " << theta.n_cols << " function evaluations." << std::endl;

  return outputValues;
}

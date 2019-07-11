#include <RcppEnsmallen.h>
#include "besselLogLikelihood.h"

// [[Rcpp::export]]
double EvaluateBessel(
    const arma::vec &p,
    const arma::mat &X,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double amplitude1 = 0,
    const double amplitude2 = 0,
    const double amplitude12 = 0,
    const double alpha1 = 0,
    const double alpha2 = 0,
    const double alpha12 = 0)
{
  // Construct the objective function.
  BesselLogLikelihood logLik;

  logLik.SetInputs(X, labels, lb, ub);

  if (alpha1 > 0)
    logLik.SetFirstAlpha(alpha1);

  if (amplitude1 > 0)
    logLik.SetFirstAmplitude(amplitude1);

  if (alpha2 > 0)
    logLik.SetSecondAlpha(alpha2);

  if (amplitude2 > 0)
    logLik.SetSecondAmplitude(amplitude2);

  if (alpha12 > 0)
    logLik.SetCrossAlpha(alpha12);

  if (amplitude12 > 0)
    logLik.SetCrossAmplitude(amplitude12);

  arma::mat params(p.n_elem, 1);
  for (unsigned int i = 0;i < p.n_elem;++i)
    params[i] = p[i];
  return logLik.Evaluate(params);
}

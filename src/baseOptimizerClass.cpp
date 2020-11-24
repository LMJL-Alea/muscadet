#include "baseOptimizerClass.h"

void BaseOptimizerFunction::TransformScaledToUnscaledParameters(arma::vec &parameters,
                                                                const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer)
{
  parameters(0) = likelihoodPointer->RetrieveAlphaFromParameters(
    parameters(0),
    likelihoodPointer->GetFirstIntensity(),
    likelihoodPointer->GetDomainDimension()
  );

  if (parameters.n_elem == 1)
    return;

  likelihoodPointer->SetFirstAlpha(parameters(0));

  parameters(1) = likelihoodPointer->RetrieveAlphaFromParameters(
    parameters(1),
    likelihoodPointer->GetSecondIntensity(),
    likelihoodPointer->GetDomainDimension()
  );

  likelihoodPointer->SetSecondAlpha(parameters(1));

  parameters(2) = likelihoodPointer->GetCrossAlphaLowerBound() / parameters(2);

  double crossAmplitude = std::sqrt(std::max(likelihoodPointer->GetSquaredCrossAmplitudeUpperBound() * parameters(3), 0.0));
  double crossIntensity = likelihoodPointer->RetrieveIntensityFromParameters(
    crossAmplitude,
    parameters(2),
    likelihoodPointer->GetDomainDimension()
  );
  parameters(3) = crossIntensity / (likelihoodPointer->GetFirstIntensity() * likelihoodPointer->GetSecondIntensity());
}

double BaseOptimizerFunction::MaximizeLikelihoodCostFunction(unsigned n,
                                                             const double *x,
                                                             double *grad,
                                                             void *data)
{
  MaximizeLikelihoodData *d = (MaximizeLikelihoodData *) data;

  arma::vec params(n);

  for (unsigned int i = 0;i < n;++i)
    params(i) = x[i];

  TransformScaledToUnscaledParameters(params, d->likelihoodPointer);

  return d->likelihoodPointer->GetValue(params);
}

double BaseOptimizerFunction::MaximizeLikelihood(arma::vec &parameters,
                                                 const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer)
{
  unsigned int numberOfParameters = likelihoodPointer->GetNumberOfParameters();
  m_LowerBounds.resize(numberOfParameters);
  m_UpperBounds.resize(numberOfParameters);
  m_LowerBounds.fill(arma::datum::eps);
  m_UpperBounds.fill(1.0 - arma::datum::eps);
  nlopt_opt optimizer = this->GetOptimizer(numberOfParameters);

  // Parameter initialization: one might want to do something smarter here such
  // as PCF initialization
  // parameters = likelihoodPointer->GetInitialPoint();
  parameters.fill(0.5);

  MaximizeLikelihoodData extraData;
  extraData.likelihoodPointer = likelihoodPointer;

  nlopt_set_lower_bounds(optimizer, &(m_LowerBounds[0]));
  nlopt_set_upper_bounds(optimizer, &(m_UpperBounds[0]));

  nlopt_set_min_objective(optimizer, this->MaximizeLikelihoodCostFunction, &extraData);
  nlopt_set_xtol_rel(optimizer, m_ParameterRelativeTolerance);

  double fVal;
  int exitCode = nlopt_optimize(optimizer, &(parameters[0]), &fVal);

  nlopt_destroy(optimizer);

  if (exitCode < 0)
  {
    Rcpp::Rcout << "Function value:   " << fVal << std::endl;
    Rcpp::Rcout << "Parameter values: " << parameters << std::endl;
    Rcpp::Rcout << "Lower bounds:     " << m_LowerBounds.as_row() << std::endl;
    Rcpp::Rcout << "Upper bounds:     " << m_UpperBounds.as_row() << std::endl;
    Rcpp::stop("NLOPT optimization failed.");
  }

  this->TransformScaledToUnscaledParameters(parameters, likelihoodPointer);

  return fVal;
}

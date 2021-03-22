#include "baseOptimizerClass.h"

void BaseOptimizerFunction::TransformScaledToUnscaledParameters(arma::vec &parameters,
                                                                const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer)
{
  if (parameters.n_elem == 2)
  {
    // Retrieve k12max and k12min
    double k12max = likelihoodPointer->GetCrossAmplitudeUpperBound();
    double k12min = parameters(1) * k12max;

    // Retrieve k12
    if (k12min < likelihoodPointer->m_ZeroValue || k12max - k12min < likelihoodPointer->m_ZeroValue)
      parameters(0) = k12min;
    else
    {
      parameters(0) *= (k12max - k12min);
      parameters(0) += k12min;
    }

    // retrieve tau
    double rho1 = likelihoodPointer->GetFirstIntensity();
    double rho2 = likelihoodPointer->GetSecondIntensity();
    unsigned int dimension = likelihoodPointer->GetDomainDimension();
    parameters(1) = k12min / likelihoodPointer->RetrieveAmplitudeFromParameters(std::sqrt(rho1 * rho2), likelihoodPointer->GetCrossAlphaLowerBound(), dimension);

    return;
  }

  parameters(0) = likelihoodPointer->RetrieveAlphaFromParameters(
    parameters(0),
    likelihoodPointer->GetFirstIntensity(),
    likelihoodPointer->GetDomainDimension()
  );
  likelihoodPointer->SetFirstAlpha(parameters(0));

  if (parameters.n_elem == 1)
    return;

  parameters(1) = likelihoodPointer->RetrieveAlphaFromParameters(
    parameters(1),
    likelihoodPointer->GetSecondIntensity(),
    likelihoodPointer->GetDomainDimension()
  );
  likelihoodPointer->SetSecondAlpha(parameters(1));

  // set crossing parameters

  // Retrieve k12max and k12min
  double k12max = likelihoodPointer->GetCrossAmplitudeUpperBound();
  double k12min = parameters(3) * k12max;

  // Retrieve k12
  if (k12min < likelihoodPointer->m_ZeroValue || k12max - k12min < likelihoodPointer->m_ZeroValue)
    parameters(2) = k12min;
  else
  {
    parameters(2) *= (k12max - k12min);
    parameters(2) += k12min;
  }

  // retrieve tau
  double rho1 = likelihoodPointer->GetFirstIntensity();
  double rho2 = likelihoodPointer->GetSecondIntensity();
  unsigned int dimension = likelihoodPointer->GetDomainDimension();
  parameters(3) = k12min / likelihoodPointer->RetrieveAmplitudeFromParameters(std::sqrt(rho1 * rho2), likelihoodPointer->GetCrossAlphaLowerBound(), dimension);
}

void BaseOptimizerFunction::TransformUnscaledToScaledParameters(arma::vec &parameters,
                                                                const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer)
{
  if (parameters.n_elem == 2)
  {
    // Scale k12max and k12min
    double k12max = likelihoodPointer->GetCrossAmplitudeUpperBound();
    double rho1 = likelihoodPointer->GetFirstIntensity();
    double rho2 = likelihoodPointer->GetSecondIntensity();
    unsigned int dimension = likelihoodPointer->GetDomainDimension();
    double k12min = likelihoodPointer->RetrieveAmplitudeFromParameters(parameters(1) * std::sqrt(rho1 * rho2), likelihoodPointer->GetCrossAlphaLowerBound(), dimension);

    if (parameters(0) < k12min)
      parameters(0) = k12min;
    if (parameters(0) > k12max)
      parameters(0) = k12max;

    // Scale k12
    if (k12min < likelihoodPointer->m_ZeroValue || k12max - k12min < likelihoodPointer->m_ZeroValue)
      parameters(0) = 0.0;
    else
    {
      parameters(0) -= k12min;
      parameters(0) /= (k12max - k12min);
    }

    // Scale tau
    parameters(1) = k12min / k12max;

    return;
  }

  likelihoodPointer->SetFirstAlpha(parameters(0));
  parameters(0) = likelihoodPointer->RetrieveAmplitudeFromParameters(
    likelihoodPointer->GetFirstIntensity(),
    parameters(0),
    likelihoodPointer->GetDomainDimension()
  );

  if (parameters.n_elem == 1)
    return;

  likelihoodPointer->SetSecondAlpha(parameters(1));
  parameters(1) = likelihoodPointer->RetrieveAmplitudeFromParameters(
    likelihoodPointer->GetSecondIntensity(),
    parameters(1),
    likelihoodPointer->GetDomainDimension()
  );

  // Scale k12max and k12min
  double k12max = likelihoodPointer->GetCrossAmplitudeUpperBound();
  double rho1 = likelihoodPointer->GetFirstIntensity();
  double rho2 = likelihoodPointer->GetSecondIntensity();
  unsigned int dimension = likelihoodPointer->GetDomainDimension();
  double k12min = likelihoodPointer->RetrieveAmplitudeFromParameters(parameters(3) * std::sqrt(rho1 * rho2), likelihoodPointer->GetCrossAlphaLowerBound(), dimension);

  if (parameters(2) < k12min)
    parameters(2) = k12min;
  if (parameters(2) > k12max)
    parameters(2) = k12max;

  // Scale k12
  if (k12min < likelihoodPointer->m_ZeroValue || k12max - k12min < likelihoodPointer->m_ZeroValue)
    parameters(2) = 0.0;
  else
  {
    parameters(2) -= k12min;
    parameters(2) /= (k12max - k12min);
  }

  // Scale tau
  parameters(3) = k12min / k12max;
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
  m_LowerBounds.fill(likelihoodPointer->m_ZeroValue);
  m_UpperBounds.fill(1.0 - likelihoodPointer->m_ZeroValue);

  if (numberOfParameters > 1)
  {
    m_LowerBounds(numberOfParameters - 2) = 0.0;
    m_UpperBounds(numberOfParameters - 2) = 1.0;
    m_LowerBounds(numberOfParameters - 1) = 0.0;
    m_UpperBounds(numberOfParameters - 1) = 1.0;
  }

  nlopt_opt optimizer = this->GetOptimizer(numberOfParameters);

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
    Rcpp::Rcout << "Exit code:   " << exitCode << std::endl;
    Rcpp::Rcout << "Function value:   " << fVal << std::endl;
    Rcpp::Rcout << "Parameter values: " << parameters.as_row() << std::endl;
    Rcpp::Rcout << (parameters[1] >= m_LowerBounds[1]) << " " << (parameters[1] <= m_UpperBounds[1]) << std::endl;
    Rcpp::Rcout << "Lower bounds:     " << m_LowerBounds.as_row() << std::endl;
    Rcpp::Rcout << "Upper bounds:     " << m_UpperBounds.as_row() << std::endl;
    Rcpp::stop("NLOPT optimization failed.");
  }

  this->TransformScaledToUnscaledParameters(parameters, likelihoodPointer);

  return fVal;
}

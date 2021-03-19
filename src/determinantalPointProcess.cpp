#include "determinantalPointProcess.h"
#include "sharedFactoryClass.h"
#include "gaussLogLikelihood.h"
#include "besselLogLikelihood.h"
#include "bobyqaOptimizerClass.h"
#include "cobylaOptimizerClass.h"
#include "neldermeadOptimizerClass.h"

Rcpp::NumericVector DeterminantalPointProcess::FormatVectorForOutput(const arma::vec &parameters) const
{
  arma::vec inputVector;

  if (parameters.n_elem == 1)
  {
    inputVector.set_size(2);
    inputVector(0) = m_LikelihoodPointer->GetFirstIntensity();
    inputVector(1) = parameters(0);
  }
  else
  {
    inputVector.set_size(7);
    inputVector(0) = m_LikelihoodPointer->GetFirstIntensity();
    inputVector(1) = m_LikelihoodPointer->GetSecondIntensity();
    inputVector(2) = (parameters.n_elem == 2) ? m_LikelihoodPointer->GetFirstAlpha() : parameters(0);
    inputVector(3) = (parameters.n_elem == 2) ? m_LikelihoodPointer->GetSecondAlpha() : parameters(1);
    inputVector(4) = (parameters.n_elem == 2) ? parameters(0) : parameters(2);
    inputVector(5) = (parameters.n_elem == 2) ? parameters(1) : parameters(3);

    if (inputVector(5) < m_LikelihoodPointer->m_ZeroValue)
      inputVector(6) = NA_REAL;
    else
    {
      double rho12Value = m_LikelihoodPointer->GetCrossIntensity();
      unsigned int dimValue = m_LikelihoodPointer->GetDomainDimension();
      inputVector(6) = m_LikelihoodPointer->RetrieveAlphaFromParameters(inputVector(4), rho12Value, dimValue);
    }
  }

  Rcpp::NumericVector outputVector = Rcpp::wrap(inputVector);
  if (outputVector.size() == 2)
    outputVector.names() = Rcpp::CharacterVector({"rho", "alpha"});
  else
    outputVector.names() = Rcpp::CharacterVector({"rho1", "rho2", "alpha1", "alpha2", "k12", "tau", "alpha12"});
  outputVector.attr("dim") = R_NilValue;
  return outputVector;
}

void DeterminantalPointProcess::SetLikelihoodModel(const std::string &val)
{
  SharedFactory<BaseLogLikelihood> likelihoodFactory;
  likelihoodFactory.Register<GaussLogLikelihood>("gauss");
  likelihoodFactory.Register<BesselLogLikelihood>("bessel");

  m_LikelihoodPointer = likelihoodFactory.Instantiate(val);

  if (!m_LikelihoodPointer)
    Rcpp::stop("The required model is not available.");
}

void DeterminantalPointProcess::SetOptimizer(const std::string &val)
{
  SharedFactory<BaseOptimizerFunction> optimizerFactory;
  optimizerFactory.Register<BobyqaOptimizerFunction>("bobyqa");
  optimizerFactory.Register<CobylaOptimizerFunction>("cobyla");
  optimizerFactory.Register<NeldermeadOptimizerFunction>("neldermead");

  m_OptimizerPointer = optimizerFactory.Instantiate(val);

  if (!m_OptimizerPointer)
    Rcpp::stop("The required optimizer is not available.");
}

Rcpp::List DeterminantalPointProcess::Fit(const arma::mat &points,
                                          const arma::vec &lower_bound,
                                          const arma::vec &upper_bound,
                                          const Rcpp::DataFrame &nd_grid,
                                          const Rcpp::Nullable<arma::vec> &init,
                                          const Rcpp::Nullable<arma::uvec> &marks,
                                          const Rcpp::Nullable<arma::vec> &marginal_parameters,
                                          const unsigned int num_threads,
                                          const unsigned int N,
                                          const unsigned int verbose_level) const
{
  // Setup data in likelihood
  arma::uvec pointMarks;

  if (marks.isNull())
  {
    pointMarks.set_size(points.n_rows);
    pointMarks.fill(1);
  }
  else
    pointMarks = Rcpp::as<arma::uvec>(marks);

  m_LikelihoodPointer->SetInputData(points, lower_bound, upper_bound, pointMarks, nd_grid, N, marginal_parameters);
  m_LikelihoodPointer->SetNumberOfThreads(num_threads);
  m_LikelihoodPointer->SetVerboseLevel(verbose_level);

  arma::vec parameters;

  if (init.isNull())
  {
    parameters.set_size(m_LikelihoodPointer->GetNumberOfParameters());
    parameters.fill(0.5);
  }
  else
  {
    parameters = Rcpp::as<arma::vec>(init);
    m_OptimizerPointer->TransformUnscaledToScaledParameters(parameters, m_LikelihoodPointer);
  }

  double minValue = m_OptimizerPointer->MaximizeLikelihood(parameters, m_LikelihoodPointer);

  return Rcpp::List::create(
    Rcpp::Named("par") = this->FormatVectorForOutput(parameters),
    Rcpp::Named("value") = minValue
  );
}

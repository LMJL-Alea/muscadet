#include "determinantalPointProcess.h"
#include "sharedFactoryClass.h"
#include "gaussLogLikelihood.h"
#include "besselLogLikelihood.h"
#include "bobyqaOptimizerClass.h"
#include "cobylaOptimizerClass.h"
#include "neldermeadOptimizerClass.h"

Rcpp::NumericVector DeterminantalPointProcess::FormatVectorForOutput(const arma::vec &inputVector) const
{
  Rcpp::NumericVector outputVector = Rcpp::wrap(inputVector);
  if (outputVector.size() == 1)
    outputVector.names() = Rcpp::CharacterVector({"alpha"});
  else
    outputVector.names() = Rcpp::CharacterVector({"alpha1", "alpha2", "alpha12", "tau"});
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
    bool validInit = true;
    for (unsigned int i = 0;i < parameters.n_elem;++i)
    {
      if (R_IsNA(parameters[i]))
      {
        validInit = false;
        break;
      }
    }
    if (!validInit)
      parameters.fill(0.5);
  }

  double minValue = m_OptimizerPointer->MaximizeLikelihood(parameters, m_LikelihoodPointer);

  if (marginal_parameters.isNotNull())
  {
    arma::vec workParams = parameters;
    parameters.set_size(4);
    parameters(0) = m_LikelihoodPointer->GetFirstAlpha();
    parameters(1) = m_LikelihoodPointer->GetSecondAlpha();
    parameters(2) = workParams(0);
    parameters(3) = workParams(1);
  }

  return Rcpp::List::create(
    Rcpp::Named("par") = this->FormatVectorForOutput(parameters),
    Rcpp::Named("value") = minValue
  );
}

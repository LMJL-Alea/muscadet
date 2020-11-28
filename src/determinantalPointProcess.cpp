#include "determinantalPointProcess.h"
#include "sharedFactoryClass.h"
#include "gaussLogLikelihood.h"
#include "besselLogLikelihood.h"
#include "bobyqaOptimizerClass.h"

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
    Rcpp::stop("The requred model is not available.");
}

void DeterminantalPointProcess::SetOptimizer(const std::string &val)
{
  SharedFactory<BaseOptimizerFunction> optimizerFactory;
  optimizerFactory.Register<BobyqaOptimizerFunction>("bobyqa");

  m_OptimizerPointer = optimizerFactory.Instantiate(val);

  if (!m_OptimizerPointer)
    Rcpp::stop("The requred optimizer is not available.");
}

Rcpp::List DeterminantalPointProcess::Fit(const arma::mat &points,
                                          const arma::vec &lower_bound,
                                          const arma::vec &upper_bound,
                                          const Rcpp::DataFrame &nd_grid,
                                          const Rcpp::Nullable<arma::uvec> &marks,
                                          const unsigned int num_threads,
                                          const unsigned int N,
                                          const bool use_verbose) const
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

  m_LikelihoodPointer->SetInputData(points, lower_bound, upper_bound, pointMarks, nd_grid, N);
  m_LikelihoodPointer->SetNumberOfThreads(num_threads);
  m_LikelihoodPointer->SetUseVerbose(use_verbose);

  arma::vec parameters(m_LikelihoodPointer->GetNumberOfParameters());
  double minValue = m_OptimizerPointer->MaximizeLikelihood(parameters, m_LikelihoodPointer);

  return Rcpp::List::create(
    Rcpp::Named("par") = this->FormatVectorForOutput(parameters),
    Rcpp::Named("value") = minValue
  );
}

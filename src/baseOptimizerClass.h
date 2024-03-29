#ifndef BASEOPTIMIZERCLASS_H
#define BASEOPTIMIZERCLASS_H

#include "baseLogLikelihood.h"

#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <memory>

class BaseOptimizerFunction
{
public:
  struct MaximizeLikelihoodData
  {
    std::shared_ptr<BaseLogLikelihood> likelihoodPointer;
  };

  BaseOptimizerFunction()
  {
    m_ParameterRelativeTolerance = std::numeric_limits<double>::epsilon();
    m_MaximumNumberOfFunctionEvaluations = 1000000;
    m_LowerBounds.reset();
    m_UpperBounds.reset();
  }

  virtual ~BaseOptimizerFunction() {}

  virtual nlopt_opt GetOptimizer(const unsigned int numberOfParameters) = 0;

  double MaximizeLikelihood(
      arma::vec &parameters,
      const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer
  );

  static void TransformScaledToUnscaledParameters(
      arma::vec &parameters,
      const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer
  );
  static void TransformUnscaledToScaledParameters(
      arma::vec &parameters,
      const std::shared_ptr<BaseLogLikelihood> &likelihoodPointer
  );

protected:
  static double MaximizeLikelihoodCostFunction(
      unsigned n,
      const double *x,
      double *grad,
      void *data
  );

private:
  double m_ParameterRelativeTolerance;
  unsigned int m_MaximumNumberOfFunctionEvaluations;
  arma::vec m_LowerBounds;
  arma::vec m_UpperBounds;
};

#endif /* BASEOPTIMIZERCLASS_H */

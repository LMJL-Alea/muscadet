#pragma once

#include "baseLogLikelihood.h"

class GaussianLogLikelihood : public BaseLogLikelihood
{
protected:
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters();
  double GetIntegral();
  double GetLogDeterminant();
};

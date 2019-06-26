#pragma once

#include "baseLogLikelihood.h"

class GaussianLogLikelihood : public BaseLogLikelihood
{
protected:
  bool CheckModelParameters();
  double GetIntegral();
  double GetLogDeterminant();
};

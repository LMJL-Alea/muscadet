#pragma once

#include "baseLogLikelihood.h"

class GaussLogLikelihood : public BaseLogLikelihood
{
public:
  void SetFirstAlpha(const double x);
  void SetSecondAlpha(const double x);
  void SetCrossAlpha(const double x);
  void SetFirstIntensity(const double x);
  void SetSecondIntensity(const double x);
  void SetCrossIntensity(const double x);

private:
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters(const arma::mat &params);
  double EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha);
  double GetIntegral();
};

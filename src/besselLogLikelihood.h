#pragma once

#include "baseLogLikelihood.h"

class BesselLogLikelihood : public BaseLogLikelihood
{
public:
  void SetFirstAlpha(const double x);
  void SetSecondAlpha(const double x);
  void SetCrossAlpha(const double x);
  void SetFirstAmplitude(const double x);
  void SetSecondAmplitude(const double x);
  void SetCrossAmplitude(const double x);

private:
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters(const arma::mat &params);
  double EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha);
  double GetIntegral();
};

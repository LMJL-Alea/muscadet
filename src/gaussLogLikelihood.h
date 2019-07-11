#include "baseLogLikelihood.h"

class GaussLogLikelihood : public BaseLogLikelihood
{
private:
  double EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha);
  double GetIntegral();
  double RetrieveIntensityFromParameters(const double amplitude, const double alpha);
  bool EvaluateAlphaConstraint();
};

#include "baseLogLikelihood.h"

class GaussLogLikelihood : public BaseLogLikelihood
{
private:
  double EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha, const unsigned int dimension);
  double RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension);
  double RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension);
  bool EvaluateAlphaConstraint(const double firstAlpha, const double secondAlpha, const double crossAlpha);
  static double GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension);
  KFunctionType GetKFunction();
};

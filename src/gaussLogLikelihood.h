#include "baseLogLikelihood.h"

class GaussLogLikelihood : public BaseLogLikelihood
{
public:
  double RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension);
  double RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension);
  double RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension);

private:
  double EvaluateLFunction(
      const double sqDist,
      const double amplitude,
      const double amplitude12,
      const double alpha,
      const double alpha12,
      const double l12Value,
      const unsigned int dimension
  );
  double EvaluateL12Function(
      const double sqDist,
      const double amplitude1,
      const double amplitude2,
      const double amplitude12,
      const double alpha12,
      const unsigned int dimension
  );
  bool EvaluateAlphaConstraint(const double firstAlpha, const double secondAlpha, const double crossAlpha);
  static double GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension);
  KFunctionType GetKFunction();
};

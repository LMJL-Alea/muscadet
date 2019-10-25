#include "baseLogLikelihood.h"

class BesselLogLikelihood : public BaseLogLikelihood
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
      const double l12Value,
      const unsigned int dimension
  );
  double EvaluateL12Function(
      const double sqDist,
      const double amplitude1,
      const double amplitude2,
      const double amplitude12,
      const double alpha12inv,
      const unsigned int dimension
  );
  double GetCrossAlphaUpperBound();
  static double GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension, const bool cross = false);
  KFunctionType GetKFunction();
};

#ifndef BESSELLOGLIKELIHOOD_H
#define BESSELLOGLIKELIHOOD_H

#include "baseLogLikelihood.h"

class BesselLogLikelihood : public BaseLogLikelihood
{
public:
  double RetrieveIntensityFromParameters(
      const double amplitude,
      const double alpha,
      const unsigned int dimension
  );
  double RetrieveAlphaFromParameters(
      const double amplitude,
      const double intensity,
      const unsigned int dimension
  );
  double RetrieveAmplitudeFromParameters(
      const double intensity,
      const double alpha,
      const unsigned int dimension
  );

private:
  double GetCrossAlphaLowerBound() const;
  double GetK11Value(const double squaredNorm);
  double GetK12Value(const double squaredNorm);
  double GetK22Value(const double squaredNorm);
};

#endif /* BESSELLOGLIKELIHOOD_H */

#include "besselLogLikelihood.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

double BesselLogLikelihood::GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension)
{
  double workValue = 2.0 * M_PI *M_PI * radius * radius;
  bool inSupport = (alpha * alpha * workValue < dimension);
  return (inSupport) ? amplitude : 0.0;
}

BesselLogLikelihood::KFunctionType BesselLogLikelihood::GetKFunction()
{
  return this->GetFourierKernel;
}

bool BesselLogLikelihood::EvaluateAlphaConstraint(const double firstAlpha, const double secondAlpha, const double crossAlpha)
{
  return crossAlpha < std::max(firstAlpha, secondAlpha);
}

double BesselLogLikelihood::EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha, const unsigned int dimension)
{
  double tmpVal = std::sqrt(2.0 * dimension * sqDist) / alpha;

  if (tmpVal < std::sqrt(std::numeric_limits<double>::epsilon()))
    return intensity / (1.0 - amplitude);

  double resVal = amplitude / (1.0 - amplitude);
  resVal *= std::pow((double)dimension / (tmpVal * M_PI * alpha * alpha), (double)dimension / 2.0);
  resVal *= boost::math::cyl_bessel_j((double)dimension / 2.0, tmpVal);
  return resVal;
}

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  double inPowerValue = 2.0 * M_PI * alpha * alpha / dimension;
  double powerValue = (double)dimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  return amplitude / denomValue;
}

double BesselLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  return std::pow(amplitude / (intensity * boost::math::tgamma(1.0 + (double)dimension / 2.0)), 1.0 / (double)dimension) * std::sqrt((double)dimension / (2.0 * M_PI));
}

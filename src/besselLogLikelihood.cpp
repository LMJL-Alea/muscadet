#include "besselLogLikelihood.h"
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

double BesselLogLikelihood::EvaluateLFunction(
    const double sqDist,
    const double amplitude,
    const double amplitude12,
    const double alpha,
    const double l12Value,
    const unsigned int dimension)
{
  double denomValue = 1.0 - amplitude;
  double firstTerm = std::pow((double)dimension / (2.0 * M_PI * alpha * alpha), (double)dimension / 2.0);
  double secondTerm = this->GetBesselJRatio(sqDist, alpha, dimension);
  return (amplitude12 * l12Value + amplitude * firstTerm * secondTerm) / denomValue;
}

double BesselLogLikelihood::EvaluateL12Function(
    const double sqDist,
    const double amplitude1,
    const double amplitude2,
    const double amplitude12,
    const double alpha12,
    const unsigned int dimension)
{
  double firstTerm = std::pow((double)dimension / (2.0 * M_PI * alpha12 * alpha12), (double)dimension / 2.0);
  double secondTerm = amplitude12 / ((1.0 - amplitude1) * (1.0 - amplitude2) - amplitude12 * amplitude12);
  double thirdTerm = this->GetBesselJRatio(sqDist, alpha12, dimension);
  return firstTerm * secondTerm * thirdTerm;
}

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = M_PI * alpha * alpha / order;
  double denomValue = std::pow(inPowerValue, order) * boost::math::tgamma(1.0 + order);
  return amplitude / denomValue;
}

double BesselLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  if (intensity < std::sqrt(std::numeric_limits<double>::epsilon()))
    return std::sqrt(std::numeric_limits<double>::epsilon());
  double order = (double)dimension / 2.0;
  return std::pow(amplitude / (intensity * boost::math::tgamma(1.0 + order)), 2.0 * order) * std::sqrt(order / M_PI);
}

double BesselLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(M_PI * alpha * alpha / order, order) * boost::math::tgamma(1.0 + order);
}

#include "besselLogLikelihood.h"
#include <boost/math/special_functions/gamma.hpp>

double BesselLogLikelihood::GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension, const bool cross)
{
  // if cross is true, alpha is its inverse
  double workValue = std::sqrt(2.0 / (double)dimension) * M_PI * radius;
  bool inSupport = (cross) ? (workValue < alpha) : (alpha * workValue < 1.0);
  return (inSupport) ? amplitude : 0.0;
}

BesselLogLikelihood::KFunctionType BesselLogLikelihood::GetKFunction()
{
  return this->GetFourierKernel;
}

double BesselLogLikelihood::GetCrossAlphaUpperBound()
{
  return std::max(this->GetFirstAlpha(), this->GetSecondAlpha());
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
    const double alpha12inv,
    const unsigned int dimension)
{
  double firstTerm = std::pow((double)dimension * alpha12inv * alpha12inv / (2.0 * M_PI), (double)dimension / 2.0);
  double secondTerm = amplitude12 / ((1.0 - amplitude1) * (1.0 - amplitude2) - amplitude12 * amplitude12);
  double thirdTerm = this->GetBesselJRatio(sqDist, alpha12inv, dimension, true);
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
  double lim = 0.001; // std::sqrt(std::numeric_limits<double>::epsilon())
  if (intensity < lim)
    return std::sqrt(std::numeric_limits<double>::epsilon());
  double order = (double)dimension / 2.0;
  double gamma = boost::math::tgamma(1.0 + order);
  return std::pow(amplitude / (intensity * gamma), 1.0 / (2.0 * order)) * std::sqrt(order / M_PI);
}

double BesselLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(M_PI * alpha * alpha / order, order) * boost::math::tgamma(1.0 + order);
}

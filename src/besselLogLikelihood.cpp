#include "besselLogLikelihood.h"
#include <boost/math/special_functions/gamma.hpp>

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = M_PI * alpha * alpha / order;
  double denomValue = std::pow(inPowerValue, order) * boost::math::tgamma(1.0 + order);
  return amplitude / denomValue;
}

double BesselLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double gammaValue = boost::math::tgamma(1.0 + order);
  return std::pow(amplitude / (intensity * gammaValue), 1.0 / (2.0 * order)) * std::sqrt(order / M_PI);
}

double BesselLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(M_PI * alpha * alpha / order, order) * boost::math::tgamma(1.0 + order);
}

double BesselLogLikelihood::GetCrossAlphaLowerBound()
{
  return std::max(this->GetFirstAlpha(), this->GetSecondAlpha());
}

double BesselLogLikelihood::GetK11Value(const double squaredNorm)
{
  return 0.0;
}

double BesselLogLikelihood::GetK12Value(const double squaredNorm)
{
  return 0.0;
}

double BesselLogLikelihood::GetK22Value(const double squaredNorm)
{
  return 0.0;
}

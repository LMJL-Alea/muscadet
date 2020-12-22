#include "besselLogLikelihood.h"

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = M_PI * alpha * alpha / order;
  double denomValue = std::pow(inPowerValue, order) * std::tgamma(1.0 + order);
  return amplitude / denomValue;
}

double BesselLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double gammaValue = std::tgamma(1.0 + order);
  return std::pow(amplitude / (intensity * gammaValue), 1.0 / (2.0 * order)) * std::sqrt(order / M_PI);
}

double BesselLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(M_PI * alpha * alpha / order, order) * std::tgamma(1.0 + order);
}

double BesselLogLikelihood::GetSquaredCrossAlphaLowerBound() const
{
  double maxAlpha = std::max(this->GetFirstAlpha(), this->GetSecondAlpha());
  return maxAlpha * maxAlpha;
}

double BesselLogLikelihood::GetK11Value(const double squaredNorm)
{
  double alphaValue = this->GetFirstAlpha();
  if (2.0 * M_PI * M_PI * alphaValue * alphaValue * squaredNorm < (double)this->GetDomainDimension())
    return this->GetFirstAmplitude();
  return 0.0;
}

double BesselLogLikelihood::GetK12Value(const double squaredNorm)
{
  double alphaValue = this->GetCrossAlpha();
  if (2.0 * M_PI * M_PI * alphaValue * alphaValue * squaredNorm < (double)this->GetDomainDimension())
    return this->GetCrossAmplitude();
  return 0.0;
}

double BesselLogLikelihood::GetK22Value(const double squaredNorm)
{
  double alphaValue = this->GetSecondAlpha();
  if (2.0 * M_PI * M_PI * alphaValue * alphaValue * squaredNorm < (double)this->GetDomainDimension())
    return this->GetSecondAmplitude();
  return 0.0;
}

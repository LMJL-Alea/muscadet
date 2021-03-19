#include "besselLogLikelihood.h"

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = this->m_PiValue * alpha * alpha / order;
  double denomValue = std::pow(inPowerValue, order) * std::tgamma(1.0 + order);
  return amplitude / denomValue;
}

double BesselLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double gammaValue = std::tgamma(1.0 + order);
  return std::pow(amplitude / (intensity * gammaValue), 1.0 / (2.0 * order)) * std::sqrt(order / this->m_PiValue);
}

double BesselLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(this->m_PiValue * alpha * alpha / order, order) * std::tgamma(1.0 + order);
}

double BesselLogLikelihood::GetCrossAlphaLowerBound() const
{
  return std::max(this->GetFirstAlpha(), this->GetSecondAlpha());
}

double BesselLogLikelihood::GetK11Value(const double squaredNorm)
{
  double alphaValue = this->GetFirstAlpha();
  if (2.0 * this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm < (double)this->GetDomainDimension())
    return this->GetFirstAmplitude();
  return 0.0;
}

double BesselLogLikelihood::GetK12Value(const double squaredNorm)
{
  double tauValue = this->GetCorrelation();
  if (tauValue < this->m_ZeroValue)
    return 0.0;

  double kValue = this->GetCrossAmplitude();
  double rhoValue = this->GetCrossIntensity();
  unsigned int dimValue = this->GetDomainDimension();
  double alphaValue = this->RetrieveAlphaFromParameters(kValue, rhoValue, dimValue);

  if (2.0 * this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm < (double)dimValue)
    return this->GetCrossAmplitude();

  return 0.0;
}

double BesselLogLikelihood::GetK22Value(const double squaredNorm)
{
  double alphaValue = this->GetSecondAlpha();
  if (2.0 * this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm < (double)this->GetDomainDimension())
    return this->GetSecondAmplitude();
  return 0.0;
}

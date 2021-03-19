#include "gaussLogLikelihood.h"

double GaussLogLikelihood::RetrieveIntensityFromParameters(const double amplitude,
                                                           const double alpha,
                                                           const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = this->m_PiValue * alpha * alpha;
  double denomValue = std::pow(inPowerValue, order);
  return amplitude / denomValue;
}

double GaussLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  return std::pow(amplitude / intensity, 1.0 / (double)dimension) / std::sqrt(this->m_PiValue);
}

double GaussLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(this->m_PiValue * alpha * alpha, order);
}

double GaussLogLikelihood::GetCrossAlphaLowerBound() const
{
  double alpha1 = this->GetFirstAlpha();
  double alpha2 = this->GetSecondAlpha();
  return std::sqrt((alpha1 * alpha1 + alpha2 * alpha2) / 2.0);
}

double GaussLogLikelihood::GetK11Value(const double squaredNorm)
{
  double kValue = this->GetFirstAmplitude();
  double alphaValue = this->GetFirstAlpha();
  return kValue * std::exp(-this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm);
}

double GaussLogLikelihood::GetK12Value(const double squaredNorm)
{
  double tauValue = this->GetCorrelation();
  if (tauValue < this->m_ZeroValue)
    return 0.0;

  double kValue = this->GetCrossAmplitude();
  double rhoValue = this->GetCrossIntensity();
  unsigned int dimValue = this->GetDomainDimension();
  double alphaValue = this->RetrieveAlphaFromParameters(kValue, rhoValue, dimValue);

  return kValue * std::exp(-this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm);
}

double GaussLogLikelihood::GetK22Value(const double squaredNorm)
{
  double kValue = this->GetSecondAmplitude();
  double alphaValue = this->GetSecondAlpha();
  return kValue * std::exp(-this->m_PiValue * this->m_PiValue * alphaValue * alphaValue * squaredNorm);
}

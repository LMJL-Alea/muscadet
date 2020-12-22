#include "gaussLogLikelihood.h"

double GaussLogLikelihood::RetrieveIntensityFromParameters(const double amplitude,
                                                           const double alpha,
                                                           const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double inPowerValue = M_PI * alpha * alpha;
  double denomValue = std::pow(inPowerValue, order);
  return amplitude / denomValue;
}

double GaussLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  return std::pow(amplitude / intensity, 1.0 / (double)dimension) / std::sqrt(M_PI);
}

double GaussLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  return intensity * std::pow(M_PI * alpha * alpha, order);
}

double GaussLogLikelihood::GetSquaredCrossAlphaLowerBound() const
{
  double alpha1 = this->GetFirstAlpha();
  double alpha2 = this->GetSecondAlpha();
  return (alpha1 * alpha1 + alpha2 * alpha2) / 2.0;
}

double GaussLogLikelihood::GetK11Value(const double squaredNorm)
{
  double kValue = this->GetFirstAmplitude();
  double alphaValue = this->GetFirstAlpha();
  return kValue * std::exp(-M_PI * M_PI * alphaValue * alphaValue * squaredNorm);
}

double GaussLogLikelihood::GetK12Value(const double squaredNorm)
{
  double kValue = this->GetCrossAmplitude();
  double alphaValue = this->GetCrossAlpha();
  return kValue * std::exp(-M_PI * M_PI * alphaValue * alphaValue * squaredNorm);
}

double GaussLogLikelihood::GetK22Value(const double squaredNorm)
{
  double kValue = this->GetSecondAmplitude();
  double alphaValue = this->GetSecondAlpha();
  return kValue * std::exp(-M_PI * M_PI * alphaValue * alphaValue * squaredNorm);
}

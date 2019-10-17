#include "gaussLogLikelihood.h"

double GaussLogLikelihood::GetFourierKernel(const double radius, const double amplitude, const double alpha, const unsigned int dimension)
{
  double workValue = M_PI * M_PI * radius * radius;
  return amplitude * std::exp(-workValue * alpha * alpha);
}

GaussLogLikelihood::KFunctionType GaussLogLikelihood::GetKFunction()
{
  return this->GetFourierKernel;
}

bool GaussLogLikelihood::EvaluateAlphaConstraint(const double firstAlpha, const double secondAlpha, const double crossAlpha)
{
  return 2.0 * crossAlpha * crossAlpha < firstAlpha * firstAlpha + secondAlpha * secondAlpha;
}

double GaussLogLikelihood::EvaluateLFunction(
    const double sqDist,
    const double amplitude,
    const double amplitude12,
    const double alpha,
    const double l12Value,
    const unsigned int dimension)
{
  return 0.0;
}

double GaussLogLikelihood::EvaluateL12Function(
    const double sqDist,
    const double amplitude1,
    const double amplitude2,
    const double amplitude12,
    const double alpha12,
    const unsigned int dimension)
{
  return 0.0;
}

double GaussLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  return amplitude / std::pow(M_PI * alpha * alpha, (double)dimension / 2.0);
}

double GaussLogLikelihood::RetrieveAlphaFromParameters(const double amplitude, const double intensity, const unsigned int dimension)
{
  return std::pow(amplitude / intensity, 1.0 / (double)dimension) / std::sqrt(M_PI);
}

double GaussLogLikelihood::RetrieveAmplitudeFromParameters(const double intensity, const double alpha, const unsigned int dimension)
{
  return intensity * std::pow(M_PI * alpha * alpha, (double)dimension / 2.0);
}

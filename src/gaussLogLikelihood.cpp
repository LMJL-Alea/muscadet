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

double GaussLogLikelihood::EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha, const unsigned int dimension)
{
  const unsigned int N = 50;

  double sumVal = 0.0;

  for (unsigned int k = 1;k <= N;++k)
  {
    double tmpVal = std::pow((double)k, -(double)dimension / 2.0);
    tmpVal *= std::pow(amplitude, (double)k - 1.0);
    tmpVal *= std::exp(-sqDist / ((double)k * alpha * alpha));
    sumVal += tmpVal;
  }

  return intensity * sumVal;
}

double GaussLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha, const unsigned int dimension)
{
  return amplitude / std::pow(std::sqrt(M_PI) * alpha, (double)dimension);
}

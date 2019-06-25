#include "integrandFunctions.h"

double BaseIntegrand::operator()(const double &radius)
{
  arma::vec kernel = GetFourierTransformedKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);
  return (std::log1p(-lmax) + std::log1p(-lmin)) * radius;
}

void BaseIntegrand::RetrieveEigenvalues(const arma::vec &kernelMatrix, double &lambdaMax, double &lambdaMin)
{
  double k11 = kernelMatrix[0];
  double k12 = kernelMatrix[1];
  double k22 = kernelMatrix[2];
  double meanDiagonal = (k11 + k22) / 2.0;
  double addOn = std::sqrt((k11 - k22) * (k11 - k22) + std::abs(k12)) / 2.0;

  lambdaMax = meanDiagonal + addOn;
  lambdaMin = meanDiagonal - addOn;
}

arma::vec GaussianIntegrand::GetFourierTransformedKernel(const double radius)
{
  arma::vec out(3);

  double alpha1 = std::exp(m_Parameters[0]);
  double alpha12 = std::exp(m_Parameters[1]);
  double alpha2 = std::exp(m_Parameters[2]);
  double tau = m_Parameters[3];

  double workValue = -M_PI * M_PI * radius * radius;

  out[0] = m_Intensity1 * std::pow(alpha1 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * alpha1 * alpha1);
  out[1] = tau * std::pow(alpha12 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * alpha12 * alpha12);
  out[2] = m_Intensity2 * std::pow(alpha2 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * alpha2 * alpha2);

  return out;
}

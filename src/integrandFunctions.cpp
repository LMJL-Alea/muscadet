#include "integrandFunctions.h"

double BaseIntegrand::operator()(const double radius)
{
  arma::vec kernel = GetFourierKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);
  return (std::log1p(-lmax) + std::log1p(-lmin)) * radius;
}

void BaseIntegrand::SetParameters(const arma::mat &params)
{
  m_Alpha1     = std::exp(params[0]);
  m_Alpha12    = std::exp(params[1]);
  m_Alpha2     = std::exp(params[2]);
  m_Covariance = params[3];
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

arma::vec GaussianIntegrand::GetFourierKernel(const double radius)
{
  arma::vec out(7);

  double workValue = -M_PI * M_PI * radius * radius;

  // K11
  out[0] = m_Intensity1 * std::pow(m_Alpha1 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * m_Alpha1 * m_Alpha1);
  // K12
  double tmpValue = std::pow(m_Alpha12 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * m_Alpha12 * m_Alpha12);
  out[1] = m_Covariance * tmpValue;
  // K22
  out[2] = m_Intensity2 * std::pow(m_Alpha2 * std::sqrt(M_PI), (double)m_DataDimension) * std::exp(workValue * m_Alpha2 * m_Alpha2);
  // D_a1 (K11)
  out[3] = (m_DataDimension / m_Alpha1 + 2.0 * workValue * m_Alpha1) * out[0];
  // D_a12 (K12)
  out[4] = (m_DataDimension / m_Alpha12 + 2.0 * workValue * m_Alpha12) * out[1];
  // D_tau (K12)
  out[5] = tmpValue;
  // D_a2 (K22)
  out[6] = (m_DataDimension / m_Alpha2 + 2.0 * workValue * m_Alpha2) * out[2];

  return out;
}

double GaussianAlpha1Integrand::operator()(const double radius)
{
  arma::vec kernel = GetFourierKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);

  double diffValue = kernel[0] - kernel[2];
  double sqrtValue = std::sqrt(diffValue * diffValue + kernel[1] * kernel[1]);
  double resValue = -0.5 * m_Alpha1 * kernel[3] * ((1.0 + diffValue / sqrtValue) / (1.0 - lmax) + (1.0 - diffValue / sqrtValue) / (1.0 - lmin));
  return resValue * radius;
}

double GaussianAlpha12Integrand::operator()(const double radius)
{
  arma::vec kernel = GetFourierKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);

  double diffValue = kernel[0] - kernel[2];
  double sqrtValue = std::sqrt(diffValue * diffValue + kernel[1] * kernel[1]);
  double resValue = -0.5 * m_Alpha12 * kernel[1] * kernel[4] / sqrtValue * (1.0 / (1.0 - lmax) - 1.0 / (1.0 - lmin));
  return resValue * radius;
}

double GaussianAlpha2Integrand::operator()(const double radius)
{
  arma::vec kernel = GetFourierKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);

  double diffValue = kernel[0] - kernel[2];
  double sqrtValue = std::sqrt(diffValue * diffValue + kernel[1] * kernel[1]);
  double resValue = -0.5 * m_Alpha2 * kernel[5] * ((1.0 - diffValue / sqrtValue) / (1.0 - lmax) + (1.0 + diffValue / sqrtValue) / (1.0 - lmin));
  return resValue * radius;
}

double GaussianCovarianceIntegrand::operator()(const double radius)
{
  arma::vec kernel = GetFourierKernel(radius);
  double lmax, lmin;
  RetrieveEigenvalues(kernel, lmax, lmin);

  double diffValue = kernel[0] - kernel[2];
  double sqrtValue = std::sqrt(diffValue * diffValue + kernel[1] * kernel[1]);
  double resValue = -0.5 * kernel[1] * kernel[6] / sqrtValue * (1.0 / (1.0 - lmax) - 1.0 / (1.0 - lmin));
  return resValue * radius;
}

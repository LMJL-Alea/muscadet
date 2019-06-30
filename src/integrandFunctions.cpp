#include "integrandFunctions.h"

void BaseIntegrand::Update(const double radius)
{
  m_Kernel = GetFourierKernel(radius);
  this->RetrieveEigenvalues(m_Kernel);

  m_DiffValue = m_Kernel[0] - m_Kernel[2];
  m_SqrtValue = std::sqrt(m_DiffValue * m_DiffValue + m_Kernel[1] * m_Kernel[1]);
}

double BaseIntegrand::operator()(const double radius)
{
  this->Update(radius);
  return (std::log1p(-m_LambdaMax) + std::log1p(-m_LambdaMin)) * radius;
}

double BaseIntegrand::GetDerivativeWRTFirstAlpha(const double radius)
{
  this->Update(radius);

  double resValue = -0.5 * m_Kernel[3];

  if (m_SqrtValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    resValue *= (1.0 / (1.0 - m_LambdaMax) + 1.0 / (1.0 - m_LambdaMin));
  else
    resValue *= ((1.0 + m_DiffValue / m_SqrtValue) / (1.0 - m_LambdaMax) + (1.0 - m_DiffValue / m_SqrtValue) / (1.0 - m_LambdaMin));

  return resValue * radius;
}

double BaseIntegrand::GetDerivativeWRTCrossAlpha(const double radius)
{
  this->Update(radius);

  if (m_SqrtValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    return 0.0;

  double resValue = -0.5 * m_Kernel[1] * m_Kernel[4] / m_SqrtValue * (1.0 / (1.0 - m_LambdaMax) - 1.0 / (1.0 - m_LambdaMin));
  return resValue * radius;
}

double BaseIntegrand::GetDerivativeWRTSecondAlpha(const double radius)
{
  this->Update(radius);

  double resValue = -0.5 * m_Kernel[5];
  if (m_SqrtValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    resValue *= (1.0 / (1.0 - m_LambdaMax) + 1.0 / (1.0 - m_LambdaMin));
  else
    resValue *= ((1.0 - m_DiffValue / m_SqrtValue) / (1.0 - m_LambdaMax) + (1.0 + m_DiffValue / m_SqrtValue) / (1.0 - m_LambdaMin));

  return resValue * radius;
}

double BaseIntegrand::GetDerivativeWRTCrossIntensity(const double radius)
{
  this->Update(radius);

  if (m_SqrtValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    return 0.0;

  double resValue = -0.5 * m_Kernel[1] * m_Kernel[6] / m_SqrtValue * (1.0 / (1.0 - m_LambdaMax) - 1.0 / (1.0 - m_LambdaMin));
  return resValue * radius;
}


void BaseIntegrand::RetrieveEigenvalues(const arma::vec &kernelMatrix)
{
  double k11 = kernelMatrix[0];
  double k12 = kernelMatrix[1];
  double k22 = kernelMatrix[2];
  double meanDiagonal = (k11 + k22) / 2.0;
  double addOn = std::sqrt((k11 - k22) * (k11 - k22) + k12 * k12) / 2.0;

  m_LambdaMax = meanDiagonal + addOn;
  m_LambdaMin = meanDiagonal - addOn;
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
  out[3] = (m_DataDimension + 2.0 * workValue * m_Alpha1 * m_Alpha1) * out[0];
  // D_a12 (K12)
  out[4] = (m_DataDimension + 2.0 * workValue * m_Alpha12 * m_Alpha12) * out[1];
  // D_a2 (K22)
  out[5] = (m_DataDimension + 2.0 * workValue * m_Alpha2 * m_Alpha2) * out[2];
  // D_tau (K12)
  out[6] = tmpValue;

  return out;
}

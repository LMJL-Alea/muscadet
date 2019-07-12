#include "integrandFunctions.h"

void BaseIntegrand::Update(const double radius)
{
  m_Kernel.set_size(7);
  m_Kernel.fill(0.0);
  m_Kernel[0] = m_KFunction(radius, m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
  m_Kernel[1] = m_KFunction(radius, m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);
  m_Kernel[2] = m_KFunction(radius, m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
  this->RetrieveEigenvalues(m_Kernel);

  m_DiffValue = m_Kernel[0] - m_Kernel[2];
  m_SqrtValue = std::sqrt(m_DiffValue * m_DiffValue + 4.0 * m_Kernel[1] * m_Kernel[1]);
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
  double addOn = std::sqrt((k11 - k22) * (k11 - k22) + 4.0 * k12 * k12) / 2.0;

  m_LambdaMax = meanDiagonal + addOn;
  m_LambdaMin = meanDiagonal - addOn;
}

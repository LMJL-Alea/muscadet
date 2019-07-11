#include "integrandFunctions.h"

void BaseIntegrand::Update(const double radius)
{
  m_Kernel = GetFourierKernel(radius);
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

arma::vec GaussianIntegrand::GetFourierKernel(const double radius)
{
  arma::vec out(7);

  double workValue = -1.0 * M_PI * M_PI * radius * radius;

  // K11
  out[0] = m_FirstIntensity * std::pow(m_FirstAlpha * std::sqrt(M_PI), (double)m_DomainDimension) * std::exp(workValue * m_FirstAlpha * m_FirstAlpha);
  // K12
  double tmpValue = std::pow(m_CrossAlpha * std::sqrt(M_PI), (double)m_DomainDimension) * std::exp(workValue * m_CrossAlpha * m_CrossAlpha);
  out[1] = m_CrossIntensity * tmpValue;
  // K22
  out[2] = m_SecondIntensity * std::pow(m_SecondAlpha * std::sqrt(M_PI), (double)m_DomainDimension) * std::exp(workValue * m_SecondAlpha * m_SecondAlpha);
  // D_a1 (K11)
  out[3] = (m_DomainDimension + 2.0 * workValue * m_FirstAlpha * m_FirstAlpha) * out[0];
  // D_a12 (K12)
  out[4] = (m_DomainDimension + 2.0 * workValue * m_CrossAlpha * m_CrossAlpha) * out[1];
  // D_a2 (K22)
  out[5] = (m_DomainDimension + 2.0 * workValue * m_SecondAlpha * m_SecondAlpha) * out[2];
  // D_tau (K12)
  out[6] = tmpValue;

  return out;
}

arma::vec BesselIntegrand::GetFourierKernel(const double radius)
{
  arma::vec out(7);
  out.fill(0.0);

  double workValue = 2.0 * M_PI *M_PI * radius * radius;

  // K11
  bool inSupport = (m_FirstAlpha * m_FirstAlpha * workValue < m_DomainDimension);
  out[0] = (inSupport) ? m_FirstAmplitude : 0.0;
  // K12
  inSupport = (m_CrossAlpha * m_CrossAlpha * workValue < m_DomainDimension);
  out[1] = (inSupport) ? m_CrossAmplitude : 0.0;
  // K22
  inSupport = (m_SecondAlpha * m_SecondAlpha * workValue < m_DomainDimension);
  out[2] = (inSupport) ? m_SecondAmplitude : 0.0;

  return out;
}

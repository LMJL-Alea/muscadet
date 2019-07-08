#pragma once

#include "baseLogLikelihood.h"

class GaussianLogLikelihood : public BaseLogLikelihood
{
public:
  GaussianLogLikelihood() : BaseLogLikelihood()
  {
    m_FirstAlpha = 0.0;
    m_SecondAlpha = 0.0;
    m_CrossAlpha = 0.0;
    m_FirstIntensity = 0.0;
    m_SecondIntensity = 0.0;
    m_CrossIntensity = 0.0;
  }

protected:
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters(const arma::mat &params);
  double GetIntegral();
  double GetLogDeterminant();

private:
  unsigned int GetNumberOfParameters() {return 4;}

  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
};

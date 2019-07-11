#pragma once

#include "baseLogLikelihood.h"

class BesselLogLikelihood : public BaseLogLikelihood
{
public:
  BesselLogLikelihood() : BaseLogLikelihood()
  {
    m_FirstAlpha = 1.0;
    m_SecondAlpha = 1.0;
    m_CrossAlpha = 1.0;
    m_FirstIntensity = 0.0;
    m_SecondIntensity = 0.0;
    m_CrossIntensity = 0.0;
    m_FirstAmplitude = 0.0;
    m_SecondAmplitude = 0.0;
    m_CrossAmplitude = 0.0;
    m_EstimateFirstBValue = true;
    m_EstimateSecondBValue = true;
    m_EstimateCrossBValue = true;
    m_EstimateFirstBetaValue = true;
    m_EstimateSecondBetaValue = true;
    m_EstimateCrossBetaValue = true;
  }

  void SetFirstAlpha(const double x);
  void SetSecondAlpha(const double x);
  void SetCrossAlpha(const double x);
  void SetFirstAmplitude(const double x);
  void SetSecondAmplitude(const double x);
  void SetCrossAmplitude(const double x);

private:
  unsigned int GetNumberOfParameters();
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters(const arma::mat &params);
  double EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha);
  double GetIntegral();
  double GetLogDeterminant();

  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
  bool m_EstimateFirstBValue, m_EstimateCrossBValue, m_EstimateSecondBValue;
  bool m_EstimateFirstBetaValue, m_EstimateCrossBetaValue, m_EstimateSecondBetaValue;
};

#pragma once

#include <RcppEnsmallen.h>

class BaseIntegrand
{
public:
  void SetFirstAlpha(const double x) {m_FirstAlpha = x;}
  void SetCrossAlpha(const double x) {m_CrossAlpha = x;}
  void SetSecondAlpha(const double x) {m_SecondAlpha = x;}
  void SetFirstIntensity(const double x) {m_FirstIntensity = x;}
  void SetCrossIntensity(const double x) {m_CrossIntensity = x;}
  void SetSecondIntensity(const double x) {m_SecondIntensity = x;}
  void SetDataDimension(const unsigned int d) {m_DomainDimension = d;}

  void Update(const double radius);
  double operator()(const double radius);
  double GetDerivativeWRTFirstAlpha(const double radius);
  double GetDerivativeWRTCrossAlpha(const double radius);
  double GetDerivativeWRTSecondAlpha(const double radius);
  double GetDerivativeWRTCrossIntensity(const double radius);

protected:
  virtual arma::vec GetFourierKernel(const double radius) = 0;
  void RetrieveEigenvalues(const arma::vec &kernelMatrix);

  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  unsigned int m_DomainDimension;

private:
  arma::vec m_Kernel;
  double m_LambdaMax, m_LambdaMin;
  double m_DiffValue, m_SqrtValue;
};

class GaussianIntegrand : public BaseIntegrand
{
protected:
  arma::vec GetFourierKernel(const double radius);
};

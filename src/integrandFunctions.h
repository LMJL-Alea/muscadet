#pragma once

#include <RcppEnsmallen.h>

class BaseIntegrand
{
public:
  typedef std::function<double(const double, const double, const double, const unsigned int)> KFunctionType;

  BaseIntegrand() {}
  ~BaseIntegrand() {}

  void SetKFunction(const KFunctionType f) {m_KFunction = f;}
  void SetFirstAlpha(const double x) {m_FirstAlpha = x;}
  void SetCrossAlpha(const double x) {m_CrossAlpha = x;}
  void SetSecondAlpha(const double x) {m_SecondAlpha = x;}
  void SetFirstAmplitude(const double x) {m_FirstAmplitude = x;}
  void SetCrossAmplitude(const double x) {m_CrossAmplitude = x;}
  void SetSecondAmplitude(const double x) {m_SecondAmplitude = x;}
  void SetDomainDimension(const unsigned int d) {m_DomainDimension = d;}

  double operator()(const double radius);
  double GetDerivativeWRTFirstAlpha(const double radius);
  double GetDerivativeWRTCrossAlpha(const double radius);
  double GetDerivativeWRTSecondAlpha(const double radius);
  double GetDerivativeWRTCrossIntensity(const double radius);

private:
  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
  unsigned int m_DomainDimension;

  void RetrieveEigenvalues(const arma::vec &kernelMatrix);
  void Update(const double radius);

  KFunctionType m_KFunction;
  arma::vec m_Kernel;
  double m_LambdaMax, m_LambdaMin;
  double m_DiffValue, m_SqrtValue;
};

#pragma once

#include <RcppEnsmallen.h>

class BaseIntegrand
{
public:
  void SetAlpha1(const double x) {m_Alpha1 = x;}
  void SetAlpha12(const double x) {m_Alpha12 = x;}
  void SetAlpha2(const double x) {m_Alpha2 = x;}
  void SetCovariance(const double x) {m_Covariance = x;}
  void SetIntensity1(const double x) {m_Intensity1 = x;}
  void SetIntensity2(const double x) {m_Intensity2 = x;}
  void SetDataDimension(const unsigned int d) {m_DataDimension = d;}

  void Update(const double radius);
  double operator()(const double radius);
  double GetDerivativeWRTFirstAlpha(const double radius);
  double GetDerivativeWRTCrossAlpha(const double radius);
  double GetDerivativeWRTSecondAlpha(const double radius);
  double GetDerivativeWRTCrossIntensity(const double radius);

protected:
  virtual arma::vec GetFourierKernel(const double radius) = 0;
  void RetrieveEigenvalues(const arma::vec &kernelMatrix);

  double m_Alpha1, m_Alpha2, m_Alpha12, m_Covariance;
  double m_Intensity1, m_Intensity2;
  unsigned int m_DataDimension;

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

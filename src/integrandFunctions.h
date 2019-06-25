#include <RcppEnsmallen.h>

class BaseIntegrand
{
public:
  void SetParameters(const arma::mat &params);
  void SetIntensity1(const double x) {m_Intensity1 = x;}
  void SetIntensity2(const double x) {m_Intensity2 = x;}
  void SetDataDimension(const unsigned int d) {m_DataDimension = d;}

  virtual double operator()(const double radius);

protected:
  virtual arma::vec GetFourierKernel(const double radius) = 0;
  void RetrieveEigenvalues(const arma::vec &kernelMatrix, double &lambdaMax, double &lambdaMin);
  double m_Alpha1, m_Alpha2, m_Alpha12, m_Covariance;
  double m_Intensity1, m_Intensity2;
  unsigned int m_DataDimension;
};

class GaussianIntegrand : public BaseIntegrand
{
protected:
  arma::vec GetFourierKernel(const double radius);
};

class GaussianAlpha1Integrand : public GaussianIntegrand
{
public:
  double operator()(const double radius);
};

class GaussianAlpha12Integrand : public GaussianIntegrand
{
public:
  double operator()(const double radius);
};

class GaussianAlpha2Integrand : public GaussianIntegrand
{
public:
  double operator()(const double radius);
};

class GaussianCovarianceIntegrand : public GaussianIntegrand
{
public:
  double operator()(const double radius);
};

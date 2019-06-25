#include <RcppEnsmallen.h>

class BaseIntegrand
{
public:
  void SetParameters(const arma::mat &params) {m_Parameters = params;};
  void SetIntensity1(const double x) {m_Intensity1 = x;}
  void SetIntensity2(const double x) {m_Intensity2 = x;}
  void SetDataDimension(const unsigned int d) {m_DataDimension = d;}

  double operator()(const double &radius);
  ~BaseIntegrand() {}

protected:
  virtual arma::vec GetFourierTransformedKernel(const double radius) = 0;
  void RetrieveEigenvalues(const arma::vec &kernelMatrix, double &lambdaMax, double &lambdaMin);
  arma::mat m_Parameters;
  double m_Intensity1, m_Intensity2;
  unsigned int m_DataDimension;
};

class GaussianIntegrand : public BaseIntegrand
{
protected:
  arma::vec GetFourierTransformedKernel(const double radius);
};

#include "gaussLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

bool GaussLogLikelihood::EvaluateAlphaConstraint()
{
  return 2.0 * m_CrossAlpha * m_CrossAlpha < m_FirstAlpha * m_FirstAlpha + m_SecondAlpha * m_SecondAlpha;
}

double GaussLogLikelihood::GetIntegral()
{
  typedef boost::math::quadrature::gauss_kronrod<double, 15> QuadratureType;
  const double lBound = 0.0;
  const double uBound = std::numeric_limits<double>::infinity();

  GaussianIntegrand integrand;
  integrand.SetFirstAlpha(m_FirstAlpha);
  integrand.SetCrossAlpha(m_CrossAlpha);
  integrand.SetSecondAlpha(m_SecondAlpha);
  integrand.SetCrossIntensity(m_CrossIntensity);
  integrand.SetFirstIntensity(m_FirstIntensity);
  integrand.SetSecondIntensity(m_SecondIntensity);
  integrand.SetDomainDimension(m_DomainDimension);
  auto GetIntegrandValue =              [&integrand](const double &t){return integrand(t);};
  auto GetDerivativeWRTFirstAlpha =     [&integrand](const double &t){return integrand.GetDerivativeWRTFirstAlpha(t);};
  auto GetDerivativeWRTCrossAlpha =     [&integrand](const double &t){return integrand.GetDerivativeWRTCrossAlpha(t);};
  auto GetDerivativeWRTSecondAlpha =    [&integrand](const double &t){return integrand.GetDerivativeWRTSecondAlpha(t);};
  auto GetDerivativeWRTCrossIntensity = [&integrand](const double &t){return integrand.GetDerivativeWRTCrossIntensity(t);};

  double resVal = 2.0 * M_PI * QuadratureType::integrate(GetIntegrandValue, lBound, uBound);

  m_GradientIntegral.set_size(this->GetNumberOfParameters());
  m_GradientIntegral[0] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTFirstAlpha,     lBound, uBound);
  m_GradientIntegral[1] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossAlpha,     lBound, uBound);
  m_GradientIntegral[2] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTSecondAlpha,    lBound, uBound);
  m_GradientIntegral[3] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossIntensity, lBound, uBound);

  return resVal;
}

double GaussLogLikelihood::EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha)
{
  const unsigned int N = 50;

  double resVal = 0.0;

  for (unsigned int k = 1;k <= N;++k)
  {
    double tmpVal = std::pow((double)k, -(double)m_DomainDimension / 2.0);
    double expInValue = sqDist / ((double)k * alpha * alpha);
    tmpVal *= std::pow(amplitude, (double)k - 1.0);
    tmpVal *= std::exp(-expInValue);
    resVal += intensity * tmpVal;
  }

  return resVal;
}

double GaussLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha)
{
  return amplitude / std::pow(std::sqrt(M_PI) * alpha, (double)m_DomainDimension);
}

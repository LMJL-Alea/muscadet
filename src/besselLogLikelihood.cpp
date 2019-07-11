#include "besselLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

bool BesselLogLikelihood::EvaluateAlphaConstraint()
{
  return m_CrossAlpha < std::max(m_FirstAlpha, m_SecondAlpha);
}

double BesselLogLikelihood::GetIntegral()
{
  typedef boost::math::quadrature::gauss_kronrod<double, 15> QuadratureType;
  const double lBound = 0.0;
  const double uBound = std::numeric_limits<double>::infinity();

  BesselIntegrand integrand;
  integrand.SetFirstAlpha(m_FirstAlpha);
  integrand.SetCrossAlpha(m_CrossAlpha);
  integrand.SetSecondAlpha(m_SecondAlpha);
  integrand.SetFirstAmplitude(m_FirstAmplitude);
  integrand.SetCrossAmplitude(m_CrossAmplitude);
  integrand.SetSecondAmplitude(m_SecondAmplitude);
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

double BesselLogLikelihood::EvaluateLFunction(const double sqDist, const double intensity, const double amplitude, const double alpha)
{
  double tmpVal = std::sqrt(2.0 * m_DomainDimension * sqDist) / alpha;

  if (tmpVal < std::sqrt(std::numeric_limits<double>::epsilon()))
    return intensity / (1.0 - amplitude);

  double resVal = amplitude / (1.0 - amplitude);
  resVal *= std::pow((double)m_DomainDimension / (tmpVal * M_PI * alpha * alpha), (double)m_DomainDimension / 2.0);
  resVal *= boost::math::cyl_bessel_j((double)m_DomainDimension / 2.0, tmpVal);
  return resVal;
}

double BesselLogLikelihood::RetrieveIntensityFromParameters(const double amplitude, const double alpha)
{
  double inPowerValue = 2.0 * M_PI * alpha * alpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  return amplitude / denomValue;
}

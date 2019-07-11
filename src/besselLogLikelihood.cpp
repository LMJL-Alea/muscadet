#include "besselLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

void BesselLogLikelihood::SetFirstAlpha(const double x)
{
  m_FirstAlpha = x;
  double inPowerValue = 2.0 * M_PI * m_FirstAlpha * m_FirstAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_FirstIntensity = m_FirstAmplitude / denomValue;
  m_EstimateFirstBetaValue = false;
}

void BesselLogLikelihood::SetSecondAlpha(const double x)
{
  m_SecondAlpha = x;
  double inPowerValue = 2.0 * M_PI * m_SecondAlpha * m_SecondAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_SecondIntensity = m_SecondAmplitude / denomValue;
  m_EstimateSecondBetaValue = false;
}

void BesselLogLikelihood::SetCrossAlpha(const double x)
{
  m_CrossAlpha = x;
  double inPowerValue = 2.0 * M_PI * m_CrossAlpha * m_CrossAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_CrossIntensity = m_CrossAmplitude / denomValue;
  m_EstimateCrossBetaValue = false;
}

void BesselLogLikelihood::SetFirstAmplitude(const double x)
{
  m_FirstAmplitude = x;
  double inPowerValue = 2.0 * M_PI * m_FirstAlpha * m_FirstAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_FirstIntensity = m_FirstAmplitude / denomValue;
  m_EstimateFirstBValue = false;
}

void BesselLogLikelihood::SetSecondAmplitude(const double x)
{
  m_SecondAmplitude = x;
  double inPowerValue = 2.0 * M_PI * m_SecondAlpha * m_SecondAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_SecondIntensity = m_SecondAmplitude / denomValue;
  m_EstimateSecondBValue = false;
}

void BesselLogLikelihood::SetCrossAmplitude(const double x)
{
  m_CrossAmplitude = x;
  double inPowerValue = 2.0 * M_PI * m_CrossAlpha * m_CrossAlpha / m_DomainDimension;
  double powerValue = (double)m_DomainDimension / 2.0;
  double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
  m_CrossIntensity = m_CrossAmplitude / denomValue;
  m_EstimateCrossBValue = false;
}

void BesselLogLikelihood::SetModelParameters(const arma::mat &params)
{
  m_Modified = false;

  unsigned int pos = 0;

  if (m_EstimateFirstBetaValue)
  {
    double workScalar = params[pos];

    if (m_FirstAlpha != workScalar)
    {
      m_FirstAlpha = workScalar;
      double inPowerValue = 2.0 * M_PI * m_FirstAlpha * m_FirstAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_FirstIntensity = m_FirstAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateSecondBetaValue)
  {
    double workScalar = params[pos];

    if (m_SecondAlpha != workScalar)
    {
      m_SecondAlpha = workScalar;
      double inPowerValue = 2.0 * M_PI * m_SecondAlpha * m_SecondAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_SecondIntensity = m_SecondAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateCrossBetaValue)
  {
    double workScalar = params[pos];

    if (m_CrossAlpha != workScalar)
    {
      m_CrossAlpha = workScalar;
      double inPowerValue = 2.0 * M_PI * m_CrossAlpha * m_CrossAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_CrossIntensity = m_CrossAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateFirstBValue)
  {
    double workScalar = params[pos];

    if (m_FirstAmplitude != workScalar)
    {
      m_FirstAmplitude = workScalar;
      double inPowerValue = 2.0 * M_PI * m_FirstAlpha * m_FirstAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_FirstIntensity = m_FirstAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateSecondBValue)
  {
    double workScalar = params[pos];

    if (m_SecondAmplitude != workScalar)
    {
      m_SecondAmplitude = workScalar;
      double inPowerValue = 2.0 * M_PI * m_SecondAlpha * m_SecondAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_SecondIntensity = m_SecondAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateCrossBValue)
  {
    double workScalar = params[pos];

    if (m_CrossAmplitude != workScalar)
    {
      m_CrossAmplitude = workScalar;
      double inPowerValue = 2.0 * M_PI * m_CrossAlpha * m_CrossAlpha / m_DomainDimension;
      double powerValue = (double)m_DomainDimension / 2.0;
      double denomValue = std::pow(inPowerValue, powerValue) * boost::math::tgamma(1.0 + powerValue);
      m_CrossIntensity = m_CrossAmplitude / denomValue;
      m_Modified = true;
    }

    ++pos;
  }

  Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossIntensity << " " << m_FirstAmplitude << " " << m_SecondAmplitude << " " << m_CrossAmplitude << std::endl;
}

bool BesselLogLikelihood::CheckModelParameters(const arma::mat &params)
{
  if (m_FirstAmplitude > 1.0 - m_Epsilon)
    return false;

  if (m_SecondAmplitude > 1.0 - m_Epsilon)
    return false;

  if (m_CrossAlpha < std::max(m_FirstAlpha, m_SecondAlpha))
    return false;

  if (m_CrossAmplitude * m_CrossAmplitude > std::min(m_FirstAmplitude * m_SecondAmplitude, (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude)))
    return false;

  return true;
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

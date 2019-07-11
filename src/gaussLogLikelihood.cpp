#include "gaussLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

void GaussLogLikelihood::SetFirstAlpha(const double x)
{
  m_FirstAlpha = x;
  m_FirstAmplitude = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
  m_EstimateFirstBetaValue = false;
}

void GaussLogLikelihood::SetSecondAlpha(const double x)
{
  m_SecondAlpha = x;
  m_SecondAmplitude = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
  m_EstimateSecondBetaValue = false;
}

void GaussLogLikelihood::SetCrossAlpha(const double x)
{
  m_CrossAlpha = x;
  m_CrossAmplitude  = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
  m_EstimateCrossBetaValue = false;
}

void GaussLogLikelihood::SetFirstIntensity(const double x)
{
  m_FirstIntensity = x;
  m_FirstAmplitude  = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
  m_EstimateFirstBValue = false;
}

void GaussLogLikelihood::SetSecondIntensity(const double x)
{
  m_SecondIntensity = x;
  m_SecondAmplitude  = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
  m_EstimateSecondBValue = false;
}

void GaussLogLikelihood::SetCrossIntensity(const double x)
{
  // m_CrossIntensity = x;
  // m_CrossAmplitude  = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);

  m_CrossAmplitude = x;
  m_CrossIntensity  = m_CrossAmplitude / std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);

  m_EstimateCrossBValue = false;
}

void GaussLogLikelihood::SetModelParameters(const arma::mat &params)
{
  m_Modified = false;

  unsigned int pos = 0;

  if (m_EstimateFirstBetaValue)
  {
    double workScalar = params[pos];

    if (m_FirstAlpha != workScalar)
    {
      m_FirstAlpha = workScalar;
      if (m_EstimateFirstBValue)
        m_FirstIntensity = m_FirstAmplitude / std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
      else
        m_FirstAmplitude = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
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
      if (m_EstimateSecondBValue)
        m_SecondIntensity = m_SecondAmplitude / std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
      else
        m_SecondAmplitude = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
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
      // if (m_EstimateCrossBValue)
        m_CrossIntensity = m_CrossAmplitude / std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
      // else
        // m_CrossAmplitude = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
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
      m_FirstIntensity = m_FirstAmplitude / std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
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
      m_SecondIntensity = m_SecondAmplitude / std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
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
      m_CrossIntensity = m_CrossAmplitude / std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossIntensity << " " << m_FirstAmplitude << " " << m_SecondAmplitude << " " << m_CrossAmplitude << std::endl;
}

bool GaussLogLikelihood::CheckModelParameters(const arma::mat &params)
{
  // if (params[0] < m_Epsilon)
  //   return false;
  //
  // if (2.0 * params[0] * params[0] < m_FirstAlpha * m_FirstAlpha + m_SecondAlpha * m_SecondAlpha)
  //   return false;
  //
  // if (params[1] * params[1] > std::min(m_FirstAmplitude * m_SecondAmplitude, 4.0 * (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon))
  //   return false;

  // if (m_FirstAmplitude > 1.0 - m_Epsilon)
  //   return false;
  //
  // if (m_CrossAmplitude * m_CrossAmplitude > std::min(m_FirstAmplitude * m_SecondAmplitude, (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude)))
  //   return false;

  // for (unsigned int i = 0;i < params.n_rows;++i)
  // {
  //   if (params[i] < m_Epsilon - 1.0 || params[i] > 1.0 - m_Epsilon)
  //     return false;
  // }

  return true;
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

#include "gaussianLogLikelihoodV2.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

void GaussianLogLikelihoodV2::SetFirstAlpha(const double x)
{
  m_FirstAlpha = x;
  m_FirstAmplitude = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
  m_EstimateFirstBetaValue = false;
}

void GaussianLogLikelihoodV2::SetSecondAlpha(const double x)
{
  m_SecondAlpha = x;
  m_SecondAmplitude  = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
  m_EstimateSecondBetaValue = false;
}

void GaussianLogLikelihoodV2::SetCrossAlpha(const double x)
{
  m_CrossAlpha = x;
  m_CrossAmplitude  = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
  m_EstimateCrossBetaValue = false;
}

void GaussianLogLikelihoodV2::SetFirstIntensity(const double x)
{
  m_FirstIntensity = x;
  m_FirstAmplitude  = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
  m_EstimateFirstBValue = false;
}

void GaussianLogLikelihoodV2::SetSecondIntensity(const double x)
{
  m_SecondIntensity = x;
  m_SecondAmplitude  = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
  m_EstimateSecondBValue = false;
}

void GaussianLogLikelihoodV2::SetCrossIntensity(const double x)
{
  m_CrossIntensity = x;
  m_CrossAmplitude  = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
  m_EstimateCrossBValue = false;
}

unsigned int GaussianLogLikelihoodV2::GetNumberOfParameters()
{
  unsigned int numParams = 0;

  if (m_EstimateFirstBValue)
    ++numParams;

  if (m_EstimateSecondBValue)
    ++numParams;

  if (m_EstimateCrossBValue)
    ++numParams;

  if (m_EstimateFirstBetaValue)
    ++numParams;

  if (m_EstimateSecondBetaValue)
    ++numParams;

  if (m_EstimateCrossBetaValue)
    ++numParams;

  return numParams;
}

void GaussianLogLikelihoodV2::SetModelParameters(const arma::mat &params)
{
  // Rcpp::Rcout << params.as_row() << std::endl;

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
      if (m_EstimateCrossBValue)
        m_CrossIntensity = m_CrossAmplitude / std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
      else
        m_CrossAmplitude = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
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

    //   if (m_EstimateFirstBetaValue)
        m_FirstIntensity = m_FirstAmplitude / std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
      // else
      //   m_FirstAmplitude = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha, (double)m_DomainDimension);
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
      // if (m_EstimateSecondBetaValue)
        m_SecondIntensity = m_SecondAmplitude / std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
      // else
      //   m_SecondAmplitude = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha, (double)m_DomainDimension);
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
      // if (m_EstimateCrossBetaValue)
        m_CrossIntensity = m_CrossAmplitude / std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
      // else
      //   m_CrossAmplitude = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, (double)m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossIntensity << " " << m_FirstAmplitude << " " << m_SecondAmplitude << " " << m_CrossAmplitude << std::endl;
  // Rcpp::stop("bah");
}

bool GaussianLogLikelihoodV2::CheckModelParameters(const arma::mat &params)
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

double GaussianLogLikelihoodV2::GetIntegral()
{
  typedef boost::math::quadrature::gauss_kronrod<double, 61> QuadratureType;
  const double lBound = 0.0;
  const double uBound = std::numeric_limits<double>::infinity();

  // Rcpp::Rcout << "Out integrand amplitude: " << m_CrossAmplitude << std::endl;

  GaussianIntegrand integrand;
  integrand.SetFirstAlpha(m_FirstAlpha);
  integrand.SetCrossAlpha(m_CrossAlpha);
  integrand.SetSecondAlpha(m_SecondAlpha);
  integrand.SetCrossIntensity(m_CrossIntensity);
  integrand.SetFirstIntensity(m_FirstIntensity);
  integrand.SetSecondIntensity(m_SecondIntensity);
  integrand.SetDomainDimension(m_DomainDimension);
  auto GetIntegrandValue =              [&integrand](const double &t){return integrand(t);};
  // auto GetDerivativeWRTFirstAlpha =     [&integrand](const double &t){return integrand.GetDerivativeWRTFirstAlpha(t);};
  // auto GetDerivativeWRTCrossAlpha =     [&integrand](const double &t){return integrand.GetDerivativeWRTCrossAlpha(t);};
  // auto GetDerivativeWRTSecondAlpha =    [&integrand](const double &t){return integrand.GetDerivativeWRTSecondAlpha(t);};
  // auto GetDerivativeWRTCrossIntensity = [&integrand](const double &t){return integrand.GetDerivativeWRTCrossIntensity(t);};

  double resVal = 2.0 * M_PI * QuadratureType::integrate(GetIntegrandValue, lBound, uBound);

  m_GradientIntegral.set_size(this->GetNumberOfParameters());
  // m_GradientIntegral[0] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTFirstAlpha,     lBound, uBound);
  // m_GradientIntegral[1] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossAlpha,     lBound, uBound);
  // m_GradientIntegral[2] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTSecondAlpha,    lBound, uBound);
  // m_GradientIntegral[3] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossIntensity, lBound, uBound);

  return resVal;
}

double GaussianLogLikelihoodV2::GetLogDeterminant()
{
  const unsigned int N = 50;

  arma::mat lMatrix(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv1(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv2(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv3(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv4(m_SampleSize, m_SampleSize);
  double resVal = 0.0;
  double workValue1 = 0.0;
  double workValue2 = 0.0;
  double workValue3 = 0.0;
  double workValue4 = 0.0;
  double workSign = 0.0;

  for (unsigned int i = 0;i < m_SampleSize;++i)
  {
    for (unsigned int j = i;j < m_SampleSize;++j)
    {
      double sqDist = m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j);
      unsigned int workLabel = m_PointLabels[i] + m_PointLabels[j];
      // Rcpp::Rcout << "Label pair: " << m_PointLabels[i] << " " << m_PointLabels[j] << " " << workLabel << std::endl;

      resVal = 0.0;
      workValue1 = 0.0;
      workValue2 = 0.0;
      workValue3 = 0.0;
      workValue4 = 0.0;

      for (unsigned int k = 1;k <= N;++k)
      {
        double tmpVal = std::pow((double)k, -(double)m_DomainDimension / 2.0);

        if (workLabel == 2)
        {
          double expInValue = sqDist / ((double)k * m_FirstAlpha * m_FirstAlpha);
          // Rcpp::Rcout << expInValue << std::endl;
          tmpVal *= std::pow(m_FirstAmplitude, (double)k - 1.0);
          tmpVal *= std::exp(-expInValue);
          resVal += m_FirstIntensity * tmpVal;
          workValue1 += m_FirstIntensity * tmpVal * (m_DomainDimension * (k - 1.0) + 2.0 * expInValue);

        }
        else if (workLabel == 3)
        {
          double expInValue = sqDist / ((double)k * m_CrossAlpha * m_CrossAlpha);
          tmpVal *= std::pow(m_CrossAmplitude, (double)k - 1.0);
          tmpVal *= std::exp(-expInValue);
          resVal +=  m_CrossIntensity * tmpVal;
          workValue2 += m_CrossIntensity * tmpVal * (m_DomainDimension * (k - 1.0) + 2.0 * expInValue);
          workValue4 += tmpVal * k;
        }
        else
        {
          double expInValue = sqDist / ((double)k * m_SecondAlpha * m_SecondAlpha);
          tmpVal *= std::pow(m_SecondAmplitude, (double)k - 1.0);
          tmpVal *= std::exp(-expInValue);
          resVal += m_SecondIntensity * tmpVal;
          workValue3 += m_SecondIntensity * tmpVal * (m_DomainDimension * (k - 1.0) + 2.0 * expInValue);
        }
      }

      lMatrix(i, j) = resVal;
      lMatrixDeriv1(i, j) = workValue1;
      lMatrixDeriv2(i, j) = workValue2;
      lMatrixDeriv3(i, j) = workValue3;
      lMatrixDeriv4(i, j) = workValue4;

      if (i != j)
      {
        lMatrix(j, i) = resVal;
        lMatrixDeriv1(j, i) = workValue1;
        lMatrixDeriv2(j, i) = workValue2;
        lMatrixDeriv3(j, i) = workValue3;
        lMatrixDeriv4(j, i) = workValue4;
      }
    }
  }

  // Rcpp::stop("bite");

  arma::log_det(resVal, workSign, lMatrix);
  arma::mat lMatrixInverse = arma::inv(lMatrix);

  m_GradientLogDeterminant.set_size(this->GetNumberOfParameters());
  m_GradientLogDeterminant[0] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  m_GradientLogDeterminant[1] = arma::trace(lMatrixInverse * lMatrixDeriv2);
  m_GradientLogDeterminant[2] = arma::trace(lMatrixInverse * lMatrixDeriv3);
  m_GradientLogDeterminant[3] = arma::trace(lMatrixInverse * lMatrixDeriv4);

  return resVal;
}

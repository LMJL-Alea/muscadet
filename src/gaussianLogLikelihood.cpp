#include "gaussianLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

bool GaussianLogLikelihood::CheckModelParameters()
{
  m_ConstraintVector.set_size(this->NumConstraints());
  m_ConstraintVector.fill(0.0);
  unsigned int numViolations = 0;

  if (m_Amplitude1 >= 1.0)
  {
    m_ConstraintVector[0] = DBL_MAX;
    ++numViolations;
  }

  if (m_Amplitude2 >= 1.0)
  {
    m_ConstraintVector[1] = DBL_MAX;
    ++numViolations;
  }

  if (2.0 * m_Alpha12 * m_Alpha12 < m_Alpha1 * m_Alpha1 + m_Alpha2 * m_Alpha2)
  {
    m_ConstraintVector[2] = DBL_MAX;
    ++numViolations;
  }

  if (m_Amplitude12 * m_Amplitude12 > m_Amplitude1 * m_Amplitude2)
  {
    m_ConstraintVector[3] = DBL_MAX;
    ++numViolations;
  }

  if (m_Amplitude12 * m_Amplitude12 >= 4.0 * (1.0 - m_Amplitude1) * (1.0 - m_Amplitude2))
  {
    m_ConstraintVector[4] = DBL_MAX;
    ++numViolations;
  }

  return (numViolations == 0);
}

double GaussianLogLikelihood::GetIntegral()
{
  typedef boost::math::quadrature::gauss_kronrod<double, 15> QuadratureType;
  const double lBound = 0.0;
  const double uBound = std::numeric_limits<double>::infinity();

  GaussianIntegrand integrand;
  integrand.SetAlpha1(m_Alpha1);
  integrand.SetAlpha12(m_Alpha12);
  integrand.SetAlpha2(m_Alpha2);
  integrand.SetCovariance(m_Covariance);
  integrand.SetIntensity1(m_Intensity1);
  integrand.SetIntensity2(m_Intensity2);
  integrand.SetDataDimension(m_DataDimension);
  auto GetIntegrandValue = [&integrand](const double &t){return integrand(t);};
  auto GetDerivativeWRTFirstAlpha = [&integrand](const double &t){return integrand.GetDerivativeWRTFirstAlpha(t);};
  auto GetDerivativeWRTCrossAlpha = [&integrand](const double &t){return integrand.GetDerivativeWRTCrossAlpha(t);};
  auto GetDerivativeWRTSecondAlpha = [&integrand](const double &t){return integrand.GetDerivativeWRTSecondAlpha(t);};
  auto GetDerivativeWRTCrossIntensity = [&integrand](const double &t){return integrand.GetDerivativeWRTCrossIntensity(t);};

  double resVal = 2.0 * M_PI * QuadratureType::integrate(GetIntegrandValue, lBound, uBound);
  m_GradientIntegral[0] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTFirstAlpha, lBound, uBound);
  m_GradientIntegral[1] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossAlpha, lBound, uBound);
  m_GradientIntegral[2] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTSecondAlpha, lBound, uBound);
  m_GradientIntegral[3] = 2.0 * M_PI * QuadratureType::integrate(GetDerivativeWRTCrossIntensity, lBound, uBound);

  return resVal;
}

double GaussianLogLikelihood::GetLogDeterminant()
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
      unsigned int workLabel = m_PointLabels[i] + m_PointLabels[j];

      if (workLabel == 2)
      {
        resVal = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= N;++k)
        {
          double expInValue = m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha1 * m_Alpha1);
          double tmpVal = std::pow(m_Amplitude1, k - 1.0);
          tmpVal *= std::pow((double)k, -(double)m_DataDimension / 2.0);
          tmpVal *= std::exp(-expInValue);
          resVal += m_Intensity1 * tmpVal;
          workValue1 += m_Intensity1 * tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * expInValue);
        }
      }
      else if (workLabel == 3)
      {
        resVal = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= N;++k)
        {
          double expInValue = m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha12 * m_Alpha12);
          double tmpVal = std::pow(m_Amplitude12, k - 1.0);
          tmpVal *= std::pow((double)k, -(double)m_DataDimension / 2.0);
          tmpVal *= std::exp(-expInValue);
          resVal +=  m_Covariance * tmpVal;
          workValue2 += m_Covariance * tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * expInValue);
          workValue4 += tmpVal * k;
        }
      }
      else
      {
        resVal = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= N;++k)
        {
          double expInValue = m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha2 * m_Alpha2);
          double tmpVal = std::pow(m_Amplitude2, k - 1.0);
          tmpVal *= std::pow((double)k, -(double)m_DataDimension / 2.0);
          tmpVal *= std::exp(-expInValue);
          resVal += m_Intensity2 * tmpVal;
          workValue3 += m_Intensity2 * tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * expInValue);
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

  arma::log_det(resVal, workSign, lMatrix);
  arma::mat lMatrixInverse = arma::inv(lMatrix);

  m_GradientLogDeterminant[0] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  m_GradientLogDeterminant[1] = arma::trace(lMatrixInverse * lMatrixDeriv2);
  m_GradientLogDeterminant[2] = arma::trace(lMatrixInverse * lMatrixDeriv3);
  m_GradientLogDeterminant[3] = arma::trace(lMatrixInverse * lMatrixDeriv4);

  return resVal;
}

#include "gaussianLogLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

bool GaussianLogLikelihood::CheckModelParameters()
{
  m_ConstraintVector.set_size(this->NumConstraints());
  m_ConstraintVector.fill(0.0);

  if (m_Amplitude1 >= 1.0)
  {
    m_ConstraintVector[0] = DBL_MAX;
    return false;
  }

  if (m_Amplitude2 >= 1.0)
  {
    m_ConstraintVector[1] = DBL_MAX;
    return false;
  }

  if (2.0 * m_Alpha12 * m_Alpha12 < m_Alpha1 * m_Alpha1 + m_Alpha2 * m_Alpha2)
  {
    m_ConstraintVector[2] = DBL_MAX;
    return false;
  }

  if (m_Amplitude12 * m_Amplitude12 > m_Amplitude1 * m_Amplitude2)
  {
    m_ConstraintVector[3] = DBL_MAX;
    return false;
  }

  if (m_Amplitude12 * m_Amplitude12 >= 4.0 * (1.0 - m_Amplitude1) * (1.0 - m_Amplitude2))
  {
    m_ConstraintVector[4] = DBL_MAX;
    return false;
  }

  return true;
}

double GaussianLogLikelihood::GetIntegral()
{
  GaussianIntegrand integrand;
  integrand.SetAlpha1(m_Alpha1);
  integrand.SetAlpha12(m_Alpha12);
  integrand.SetAlpha2(m_Alpha2);
  integrand.SetCovariance(m_Covariance);
  integrand.SetIntensity1(m_Intensity1);
  integrand.SetIntensity2(m_Intensity2);
  integrand.SetDataDimension(m_DataDimension);
  auto f = [&integrand](double t){ return integrand(t); };
  double resVal = 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, 0.0, std::numeric_limits<double>::infinity());

  GaussianAlpha1Integrand integrand1;
  integrand1.SetAlpha1(m_Alpha1);
  integrand1.SetAlpha12(m_Alpha12);
  integrand1.SetAlpha2(m_Alpha2);
  integrand1.SetCovariance(m_Covariance);
  integrand1.SetIntensity1(m_Intensity1);
  integrand1.SetIntensity2(m_Intensity2);
  integrand1.SetDataDimension(m_DataDimension);
  auto f1 = [&integrand1](double t){ return integrand1(t); };
  m_GradientIntegral[0] = 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f1, 0.0, std::numeric_limits<double>::infinity());

  GaussianAlpha12Integrand integrand2;
  integrand2.SetAlpha1(m_Alpha1);
  integrand2.SetAlpha12(m_Alpha12);
  integrand2.SetAlpha2(m_Alpha2);
  integrand2.SetCovariance(m_Covariance);
  integrand2.SetIntensity1(m_Intensity1);
  integrand2.SetIntensity2(m_Intensity2);
  integrand2.SetDataDimension(m_DataDimension);
  auto f2 = [&integrand2](double t){ return integrand2(t); };
  m_GradientIntegral[1] = 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f2, 0.0, std::numeric_limits<double>::infinity());

  GaussianAlpha2Integrand integrand3;
  integrand3.SetAlpha1(m_Alpha1);
  integrand3.SetAlpha12(m_Alpha12);
  integrand3.SetAlpha2(m_Alpha2);
  integrand3.SetCovariance(m_Covariance);
  integrand3.SetIntensity1(m_Intensity1);
  integrand3.SetIntensity2(m_Intensity2);
  integrand3.SetDataDimension(m_DataDimension);
  auto f3 = [&integrand3](double t){ return integrand3(t); };
  m_GradientIntegral[2] = 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f3, 0.0, std::numeric_limits<double>::infinity());

  GaussianCovarianceIntegrand integrand4;
  integrand4.SetAlpha1(m_Alpha1);
  integrand4.SetAlpha12(m_Alpha12);
  integrand4.SetAlpha2(m_Alpha2);
  integrand4.SetCovariance(m_Covariance);
  integrand4.SetIntensity1(m_Intensity1);
  integrand4.SetIntensity2(m_Intensity2);
  integrand4.SetDataDimension(m_DataDimension);
  auto f4 = [&integrand4](double t){ return integrand4(t); };
  m_GradientIntegral[3] = 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f4, 0.0, std::numeric_limits<double>::infinity());

  return resVal;
}

double GaussianLogLikelihood::GetLogDeterminant()
{
  Rcpp::Rcout << m_Amplitude1 << std::endl;
  Rcpp::Rcout << m_Amplitude12 << std::endl;
  Rcpp::Rcout << m_Amplitude2 << std::endl;

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
        for (unsigned int k = 1;k <= 50;++k)
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
        for (unsigned int k = 1;k <= 50;++k)
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
        for (unsigned int k = 1;k <= 50;++k)
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

//' @export
// [[Rcpp::export]]
arma::mat CalcDistMat(const arma::mat &x)
{
  GaussianLogLikelihood logLik;
  logLik.SetInputs(x);
  return logLik.GetDistanceMatrix();
}

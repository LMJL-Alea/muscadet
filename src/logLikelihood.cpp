#include "logLikelihood.h"
#include "integrandFunctions.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

void BaseLogLikelihood::SetInputs(const arma::mat &points, const double volume)
{
  m_DataDimension = points.n_cols - 1;
  m_SampleSize = points.n_rows;
  m_PointLabels = points.col(m_DataDimension);
  m_DataVolume = volume;

  m_DistanceMatrix.set_size(m_SampleSize, m_SampleSize);
  m_DistanceMatrix.fill(0.0);
  m_Intensity1 = 0.0;
  m_Intensity2 = 0.0;
  arma::mat dataPoints = points.cols(0, m_DataDimension - 1);
  arma::rowvec workVec1, workVec2;
  for (unsigned int i = 0;i < m_SampleSize;++i)
  {
    if (m_PointLabels[i] == 1)
      ++m_Intensity1;
    if (m_PointLabels[i] == 2)
      ++m_Intensity2;

    workVec1 = dataPoints.row(i);

    for (unsigned int j = i + 1;j < m_SampleSize;++j)
    {
      workVec2 = dataPoints.row(j);
      double workDistance = arma::norm(workVec1 - workVec2);
      m_DistanceMatrix(i, j) = workDistance;
      m_DistanceMatrix(j, i) = workDistance;
    }
  }

  m_Intensity1 /= volume;
  m_Intensity2 /= volume;
}

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  m_Modified = false;

  double workScalar = std::exp(params[0]);
  if (m_Alpha1 != workScalar)
  {
    m_Alpha1 = workScalar;
    m_Modified = true;
  }

  workScalar = std::exp(params[1]);
  if (m_Alpha12 != workScalar)
  {
    m_Alpha12 = workScalar;
    m_Modified = true;
  }

  workScalar = std::exp(params[2]);
  if (m_Alpha2 != workScalar)
  {
    m_Alpha2 = workScalar;
    m_Modified = true;
  }

  workScalar = params[3];
  if (m_Covariance != workScalar)
  {
    m_Covariance = workScalar;
    m_Modified = true;
  }
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  double logLik = 2.0 * m_DataVolume;
  logLik += m_DataVolume * m_Integral;
  logLik += m_LogDeterminant;

  return -2.0 * logLik;
}

void BaseLogLikelihood::Gradient(const arma::mat& x, arma::mat &g)
{
  this->SetModelParameters(x);

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  unsigned int numParams = x.n_cols;
  g.set_size(numParams, 1);

  for (unsigned int i = 0;i < numParams;++i)
    g[i] = m_GradientIntegral[i] + m_GradientLogDeterminant[i];

  g *= -2.0;
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
          double tmpVal = std::pow(m_Intensity1, (double)k);
          tmpVal *= std::pow(std::sqrt(M_PI) * m_Alpha1, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha1 * m_Alpha1));
          resVal += tmpVal;
          workValue1 += tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha1 * m_Alpha1));
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
          double tmpVal = std::pow(m_Covariance, k - 1.0);
          tmpVal *= std::pow(std::sqrt(M_PI) * m_Alpha12, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha12 * m_Alpha12));
          resVal +=  m_Covariance * tmpVal;
          workValue2 += m_Covariance * tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha12 * m_Alpha12));
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
          double tmpVal = std::pow(m_Intensity2, (double)k);
          tmpVal *= std::pow(std::sqrt(M_PI) * m_Alpha2, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha2 * m_Alpha2));
          resVal += tmpVal;
          workValue3 += tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * m_Alpha2 * m_Alpha2));
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
  m_GradientLogDeterminant[1] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  m_GradientLogDeterminant[2] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  m_GradientLogDeterminant[3] = arma::trace(lMatrixInverse * lMatrixDeriv1);

  return resVal;
}

double GaussianLogLikelihood::EvaluateConstraint(const unsigned int i, const arma::mat& x)
{
  this->SetModelParameters(x);
  double tau = m_Covariance  / std::sqrt(m_Intensity1 * m_Intensity2);

  if (i == 0)
    return (m_Intensity1 * std::pow(std::sqrt(M_PI) * m_Alpha1, (double)m_DataDimension) <= 1.0) ? 0.0 : DBL_MAX;

  if (i == 1)
    return (m_Intensity2 * std::pow(std::sqrt(M_PI) * m_Alpha2, (double)m_DataDimension) <= 1.0) ? 0.0 : DBL_MAX;

  if (i == 2)
    return (2.0 * m_Alpha12 * m_Alpha12 >= m_Alpha1 * m_Alpha1 + m_Alpha2 * m_Alpha2) ? 0.0 : DBL_MAX;

  if (i == 3)
    return (tau * tau * std::pow(m_Alpha12, 2.0 * m_DataDimension) <= std::pow(m_Alpha1 * m_Alpha2, (double)m_DataDimension)) ? 0.0 : DBL_MAX;

  if (i == 4)
    return (tau * tau * std::pow(m_Alpha12, 2.0 * m_DataDimension) <= 4.0 * std::pow(m_Alpha1 * m_Alpha2, (double)m_DataDimension) * (1.0 / (m_Intensity1 * std::pow(std::sqrt(M_PI) * m_Alpha1, (double)m_DataDimension)) - 1.0) * (1.0 / (m_Intensity2 * std::pow(std::sqrt(M_PI) * m_Alpha2, (double)m_DataDimension)) - 1.0)) ? 0.0 : DBL_MAX;

  return -1.0;
}

//' @export
// [[Rcpp::export]]
arma::mat CalcDistMat(const arma::mat &x)
{
  GaussianLogLikelihood logLik;
  logLik.SetInputs(x);
  return logLik.GetDistanceMatrix();
}

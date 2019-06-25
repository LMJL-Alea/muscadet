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

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  arma::vec detGradient;
  double logLik = 2.0 * m_DataVolume;
  logLik += m_DataVolume * GetNormalizationFactor(x);
  logLik += GetLogDeterminant(x, detGradient);
  return -2.0 * logLik;
}

double GaussianLogLikelihood::GetNormalizationFactor(const arma::mat &params)
{
  GaussianIntegrand integrand;
  integrand.SetParameters(params);
  integrand.SetIntensity1(m_Intensity1);
  integrand.SetIntensity2(m_Intensity2);
  integrand.SetDataDimension(m_DataDimension);
  auto f = [&integrand](double t){ return integrand(t); };
  return 2.0 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, 0.0, std::numeric_limits<double>::infinity());
}

double GaussianLogLikelihood::GetLogDeterminant(const arma::mat &params, arma::vec &grad)
{
  // Transform unconstrained parameters into real parameters
  double alpha1 = std::exp(params[0]);
  double alpha12 = std::exp(params[1]);
  double alpha2 = std::exp(params[2]);
  double tau = params[3];

  arma::mat lMatrix(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv1(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv2(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv3(m_SampleSize, m_SampleSize);
  arma::mat lMatrixDeriv4(m_SampleSize, m_SampleSize);
  double workValue = 0.0;
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
        workValue = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= 50;++k)
        {
          double tmpVal = std::pow(m_Intensity1, (double)k);
          tmpVal *= std::pow(std::sqrt(M_PI) * alpha1, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha1 * alpha1));
          workValue += tmpVal;
          workValue1 += tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha1 * alpha1));
        }
      }
      else if (workLabel == 3)
      {
        workValue = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= 50;++k)
        {
          double tmpVal = std::pow(tau, k - 1.0);
          tmpVal *= std::pow(std::sqrt(M_PI) * alpha12, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha12 * alpha12));
          workValue +=  tau * tmpVal;
          workValue2 += tau * tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha12 * alpha12));
          workValue4 += tmpVal * k;
        }
      }
      else
      {
        workValue = 0.0;
        workValue1 = 0.0;
        workValue2 = 0.0;
        workValue3 = 0.0;
        workValue4 = 0.0;
        for (unsigned int k = 1;k <= 50;++k)
        {
          double tmpVal = std::pow(m_Intensity2, (double)k);
          tmpVal *= std::pow(std::sqrt(M_PI) * alpha2, m_DataDimension * (k - 1.0));
          tmpVal *= std::pow((double)k, -m_DataDimension / 2.0);
          tmpVal *= std::exp(-m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha2 * alpha2));
          workValue += tmpVal;
          workValue3 += tmpVal * (m_DataDimension * (m_SampleSize - 1.0) + 2.0 * m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j) / (k * alpha2 * alpha2));
        }
      }

      lMatrix(i, j) = workValue;
      lMatrixDeriv1(i, j) = workValue1;
      lMatrixDeriv2(i, j) = workValue2;
      lMatrixDeriv3(i, j) = workValue3;
      lMatrixDeriv4(i, j) = workValue4;

      if (i != j)
      {
        lMatrix(j, i) = workValue;
        lMatrixDeriv1(j, i) = workValue1;
        lMatrixDeriv2(j, i) = workValue2;
        lMatrixDeriv3(j, i) = workValue3;
        lMatrixDeriv4(j, i) = workValue4;
      }
    }
  }

  arma::log_det(workValue, workSign, lMatrix);

  grad.set_size(4);
  arma::mat lMatrixInverse = arma::inv(lMatrix);
  grad[0] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  grad[1] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  grad[2] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  grad[3] = arma::trace(lMatrixInverse * lMatrixDeriv1);

  return workValue;
}

double GaussianLogLikelihood::EvaluateConstraint(const unsigned int i, const arma::mat& x)
{
  // Transform unconstrained parameters into real parameters
  double alpha1 = std::exp(x[0]);
  double alpha12 = std::exp(x[1]);
  double alpha2 = std::exp(x[2]);
  double tau = x[3]  / std::sqrt(m_Intensity1 * m_Intensity2);

  if (i == 0)
    return (m_Intensity1 * std::pow(std::sqrt(M_PI) * alpha1, (double)m_DataDimension) <= 1.0) ? 0.0 : DBL_MAX;

  if (i == 1)
    return (m_Intensity2 * std::pow(std::sqrt(M_PI) * alpha2, (double)m_DataDimension) <= 1.0) ? 0.0 : DBL_MAX;

  if (i == 2)
    return (2.0 * alpha12 * alpha12 >= alpha1 * alpha1 + alpha2 * alpha2) ? 0.0 : DBL_MAX;

  if (i == 3)
    return (tau * tau * std::pow(alpha12, 2.0 * m_DataDimension) <= std::pow(alpha1 * alpha2, (double)m_DataDimension)) ? 0.0 : DBL_MAX;

  if (i == 4)
    return (tau * tau * std::pow(alpha12, 2.0 * m_DataDimension) <= 4.0 * std::pow(alpha1 * alpha2, (double)m_DataDimension) * (1.0 / (m_Intensity1 * std::pow(std::sqrt(M_PI) * alpha1, (double)m_DataDimension)) - 1.0) * (1.0 / (m_Intensity2 * std::pow(std::sqrt(M_PI) * alpha2, (double)m_DataDimension)) - 1.0)) ? 0.0 : DBL_MAX;

  return -1.0;
}

//' @export
// [[Rcpp::export]]
arma::mat CalcDistMat(const arma::mat &x)
{
  GaussianLogLikelihood logLik;
  logLik.SetInputs(x, 1.0);
  return logLik.GetDistanceMatrix();
}

//' @export
// [[Rcpp::export]]
double CalcInt(const arma::mat &x)
{
  GaussianLogLikelihood logLik;
  logLik.SetInputs(x, 1.0);
  arma::mat beta(4, 1, arma::fill::randn);
  return logLik.GetNormalizationFactor(beta);
}

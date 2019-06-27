#include "baseLogLikelihood.h"

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

  Rcpp::Rcout << "Intensity 1: " << m_Intensity1 << std::endl;
  Rcpp::Rcout << "Intensity 2: " << m_Intensity2 << std::endl;
  Rcpp::Rcout << "Point Dimension: " << m_DataDimension << std::endl;
  Rcpp::Rcout << "Sample size: " << m_SampleSize << std::endl;
}

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  Rcpp::Rcout << params.as_row() << std::endl;
  m_Modified = false;

  double workScalar = std::exp(params[0]);
  if (m_Alpha1 != workScalar)
  {
    m_Alpha1 = workScalar;
    m_Amplitude1  = m_Intensity1 * std::pow(std::sqrt(M_PI) * m_Alpha1,  m_DataDimension);
    m_Modified = true;
  }

  workScalar = std::exp(params[1]);
  if (m_Alpha12 != workScalar)
  {
    m_Alpha12 = workScalar;
    m_Amplitude12 = m_Covariance * std::pow(std::sqrt(M_PI) * m_Alpha12, m_DataDimension);
    m_Modified = true;
  }

  workScalar = std::exp(params[2]);
  if (m_Alpha2 != workScalar)
  {
    m_Alpha2 = workScalar;
    m_Amplitude2  = m_Intensity2 * std::pow(std::sqrt(M_PI) * m_Alpha2,  m_DataDimension);
    m_Modified = true;
  }

  workScalar = params[3];
  if (m_Covariance != workScalar)
  {
    m_Covariance = workScalar;
    m_Amplitude12 = m_Covariance * std::pow(std::sqrt(M_PI) * m_Alpha12, m_DataDimension);
    m_Modified = true;
  }
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  Rcpp::Rcout << "init evaluate" << std::endl;
  this->SetModelParameters(x);

  bool validParams = this->CheckModelParameters();
  if (!validParams)
    return DBL_MAX;

  if (m_Modified)
  {
    Rcpp::Rcout << "enter modified" << std::endl;
    m_Integral = this->GetIntegral();
    Rcpp::Rcout << "done with integral in evaluate" << std::endl;
    m_LogDeterminant = this->GetLogDeterminant();
    Rcpp::Rcout << "done with determ in evaluate" << std::endl;
  }

  Rcpp::Rcout << m_Integral << std::endl;
  Rcpp::Rcout << m_LogDeterminant << std::endl;

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
    Rcpp::stop("Non finite stuff");

  double logLik = 2.0 * m_DataVolume;
  logLik += m_DataVolume * m_Integral;
  logLik += m_LogDeterminant;

  Rcpp::Rcout << "done with evaluate" << std::endl;

  return -2.0 * logLik;
}

void BaseLogLikelihood::Gradient(const arma::mat& x, arma::mat &g)
{
  Rcpp::Rcout << "init gradient" << std::endl;
  this->SetModelParameters(x);

  unsigned int numParams = x.n_cols;
  g.set_size(numParams, 1);

  bool validParams = this->CheckModelParameters();
  if (!validParams)
  {
    g.fill(0.0);
    Rcpp::Rcout << "done with gradient" << std::endl;
    return;
  }

  if (m_Modified)
  {
    Rcpp::Rcout << "enter modified" << std::endl;
    m_Integral = this->GetIntegral();
    Rcpp::Rcout << "done with integral in gradient" << std::endl;
    m_LogDeterminant = this->GetLogDeterminant();
    Rcpp::Rcout << "done with determ in gradient" << std::endl;
  }

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
    Rcpp::stop("Non finite studd in gradient");

  Rcpp::Rcout << m_Integral << std::endl;
  Rcpp::Rcout << m_LogDeterminant << std::endl;

  for (unsigned int i = 0;i < numParams;++i)
    g[i] = m_GradientIntegral[i] + m_GradientLogDeterminant[i];

  g *= -2.0;

  Rcpp::Rcout << "done with gradient" << std::endl;
}

double BaseLogLikelihood::EvaluateConstraint(const size_t i, const arma::mat& x)
{
  Rcpp::Rcout << "init with evaluate constraint" << std::endl;
  this->SetModelParameters(x);
  this->CheckModelParameters();
  Rcpp::Rcout << "done with evaluate constraint: " << m_ConstraintVector[i] << std::endl;
  return m_ConstraintVector[i];
}

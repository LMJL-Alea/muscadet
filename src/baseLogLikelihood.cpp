#include "baseLogLikelihood.h"

void BaseLogLikelihood::SetInputs(const arma::mat &points, const double rho1, const double rho2, const double volume)
{
  m_DomainDimension = points.n_cols - 1;
  m_SampleSize = points.n_rows;
  m_PointLabels = points.col(m_DomainDimension);
  m_DomainVolume = volume;
  m_FirstIntensity = rho1;
  m_SecondIntensity = rho2;

  m_DistanceMatrix.set_size(m_SampleSize, m_SampleSize);
  m_DistanceMatrix.fill(0.0);
  arma::mat dataPoints = points.cols(0, m_DomainDimension - 1);
  arma::rowvec workVec1, workVec2;

  for (unsigned int i = 0;i < m_SampleSize;++i)
  {
    workVec1 = dataPoints.row(i);

    for (unsigned int j = i + 1;j < m_SampleSize;++j)
    {
      workVec2 = dataPoints.row(j);
      double workDistance = arma::norm(workVec1 - workVec2);
      m_DistanceMatrix(i, j) = workDistance;
      m_DistanceMatrix(j, i) = workDistance;
    }
  }

  Rcpp::Rcout << "Point Dimension: " << m_DomainDimension << std::endl;
  Rcpp::Rcout << "Sample size: " << m_SampleSize << std::endl;
}

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  m_Modified = false;

  double workScalar = std::exp(params[0]);
  if (m_FirstAlpha != workScalar)
  {
    m_FirstAlpha = workScalar;
    m_FirstAmplitude  = m_FirstIntensity * std::pow(std::sqrt(M_PI) * m_FirstAlpha,  m_DomainDimension);
    m_Modified = true;
  }

  workScalar = std::exp(params[1]);
  if (m_CrossAlpha != workScalar)
  {
    m_CrossAlpha = workScalar;
    m_CrossAmplitude = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, m_DomainDimension);
    m_Modified = true;
  }

  workScalar = std::exp(params[2]);
  if (m_SecondAlpha != workScalar)
  {
    m_SecondAlpha = workScalar;
    m_SecondAmplitude  = m_SecondIntensity * std::pow(std::sqrt(M_PI) * m_SecondAlpha,  m_DomainDimension);
    m_Modified = true;
  }

  workScalar = params[3];
  if (m_CrossIntensity != workScalar)
  {
    m_CrossIntensity = workScalar;
    m_CrossAmplitude = m_CrossIntensity * std::pow(std::sqrt(M_PI) * m_CrossAlpha, m_DomainDimension);
    m_Modified = true;
  }
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  bool validParams = this->CheckModelParameters();
  if (!validParams)
    return DBL_MAX;

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
    Rcpp::stop("Non finite stuff in evaluate");

  double logLik = 2.0 * m_DomainVolume;
  logLik += m_DomainVolume * m_Integral;
  logLik += m_LogDeterminant;

  return -2.0 * logLik;
}

void BaseLogLikelihood::Gradient(const arma::mat& x, arma::mat &g)
{
  this->SetModelParameters(x);

  unsigned int numParams = x.n_rows;
  g.set_size(numParams, 1);

  bool validParams = this->CheckModelParameters();
  if (!validParams)
  {
    g.fill(0.0);
    return;
  }

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
    Rcpp::stop("Non finite stuff in gradient");

  for (unsigned int i = 0;i < numParams;++i)
    g[i] = m_GradientIntegral[i] + m_GradientLogDeterminant[i];

  g *= -2.0;
}

double BaseLogLikelihood::EvaluateConstraint(const size_t i, const arma::mat& x)
{
  this->SetModelParameters(x);
  this->CheckModelParameters();
  return m_ConstraintVector[i];
}

void BaseLogLikelihood::GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g)
{
  g.set_size(x.n_rows, 1);
  g.fill(0.0);
}

#include "baseLogLikelihood.h"

void BaseLogLikelihood::SetNeighborhood(const unsigned int n)
{
  std::vector<int> workVec(n, 0);
  m_Neighborhood.resize(1);
  m_Neighborhood[0] = workVec;
  NeighborhoodType workList;

  for (unsigned int i = 0;i < n;++i)
  {
    workList = m_Neighborhood;
    m_Neighborhood.clear();

    for (unsigned int j = 0;j < workList.size();++j)
    {
      for (int k = -1;k <= 1;++k)
      {
        workVec = workList[j];
        workVec[i] += k;
        m_Neighborhood.push_back(workVec);
      }
    }
  }
}

void BaseLogLikelihood::SetInputs(
    const arma::mat &points,
    const arma::vec &labels,
    const arma::vec &lb,
    const arma::vec &ub,
    const double rho1,
    const double rho2)
{
  m_DomainDimension = points.n_cols;
  m_SampleSize = points.n_rows;
  m_PointLabels = labels;
  m_DomainLowerBounds = lb;
  m_DomainUpperBounds = ub;
  m_DomainVolume = 1.0;
  for (unsigned int i = 0;i < m_DomainDimension;++i)
    m_DomainVolume *= (ub[i] - lb[i]);

  m_FirstIntensity = rho1;
  m_SecondIntensity = rho2;

  this->SetNeighborhood(m_DomainDimension);
  m_DistanceMatrix.set_size(m_SampleSize, m_SampleSize);
  m_DistanceMatrix.fill(0.0);
  std::vector<arma::rowvec> trialVectors;
  arma::rowvec workVec1, workVec2;

  for (unsigned int i = 0;i < m_SampleSize;++i)
  {
    workVec1 = points.row(i);

    if (m_UsePeriodicDomain)
      trialVectors = this->GetTrialVectors(workVec1, lb, ub);

    for (unsigned int j = i + 1;j < m_SampleSize;++j)
    {
      workVec2 = points.row(j);

      double workDistance = 0.0;

      if (m_UsePeriodicDomain)
      {
        for (unsigned int k = 0;k < trialVectors.size();++k)
        {
          double testDistance = arma::norm(trialVectors[k] - workVec2);
          if (testDistance < workDistance || k == 0)
            workDistance = testDistance;
        }
      }
      else
        workDistance = arma::norm(workVec1 - workVec2);

      m_DistanceMatrix(i, j) = workDistance;
      m_DistanceMatrix(j, i) = workDistance;
    }
  }

  Rcpp::Rcout << "Point Dimension: " << m_DomainDimension << std::endl;
  Rcpp::Rcout << "Sample size: " << m_SampleSize << std::endl;
}

std::vector<arma::rowvec> BaseLogLikelihood::GetTrialVectors(const arma::rowvec &x, const arma::vec &lb, const arma::vec &ub)
{
  unsigned int numTrials = m_Neighborhood.size();
  std::vector<arma::rowvec> trialVectors(numTrials);
  std::vector<int> workNeighborhood;
  arma::rowvec workVector;

  for (unsigned int i = 0;i < numTrials;++i)
  {
    workNeighborhood = m_Neighborhood[i];
    workVector = x;
    for (unsigned int j = 0;j < m_DomainDimension;++j)
      workVector[j] += (double)workNeighborhood[j] * (ub[j] - lb[j]);
    trialVectors[i] = workVector;
  }

  return trialVectors;
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

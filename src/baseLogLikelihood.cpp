#include "baseLogLikelihood.h"

const double BaseLogLikelihood::m_Epsilon = 1.0e-4;

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

void BaseLogLikelihood::SetInputs(
    const arma::mat &points,
    const arma::uvec &labels,
    const arma::vec &lb,
    const arma::vec &ub)
{
  m_DomainDimension = points.n_cols;
  m_SampleSize = points.n_rows;
  m_PointLabels = labels;
  m_DomainVolume = 1.0;
  for (unsigned int i = 0;i < m_DomainDimension;++i)
    m_DomainVolume *= (ub[i] - lb[i]);

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

  // Rcpp::Rcout << "Domain Dimension: " << m_DomainDimension << std::endl;
  // Rcpp::Rcout << "Domain Volume: " << m_DomainVolume << std::endl;
  // Rcpp::Rcout << "Sample size: " << m_SampleSize << std::endl;
  // Rcpp::Rcout << "Point labels: " << m_PointLabels.as_row() << std::endl;
}

arma::mat BaseLogLikelihood::GetInitialPoint(const double rho1, const double rho2, const double alpha1, const double alpha2)
{
  arma::mat params(this->GetNumberOfParameters(), 1);

  // Initialize CrossAlpha to the maximum between FirstAlpha and
  // SecondAlpha which automatically satisfies constraint #3
  // double alpha12 = std::max(alpha1, alpha2);
  double alpha12 = (1.0 + m_Epsilon) * std::sqrt((alpha1 * alpha1 + alpha2 * alpha2) / 2.0);

  // Finds lower and upper bounds for CrossIntensity and
  // randomly sample within the corresponding interval
  double amp1 = rho1 * std::pow(std::sqrt(M_PI) * alpha1, (double)m_DomainDimension);
  double amp2 = rho2 * std::pow(std::sqrt(M_PI) * alpha2, (double)m_DomainDimension);
  double amp12 = std::pow(std::sqrt(M_PI) * alpha12, (double)m_DomainDimension);
  double ub = std::sqrt(amp1 * amp2) / amp12;
  ub = std::min(ub, 2.0 * std::sqrt((1.0 - amp1) * (1.0 - amp2)) / amp12 - m_Epsilon);
  Rcpp::Rcout << ub << std::endl;
  double rho12 = ub * (2.0 * arma::randu() - 1.0);

  params[0] = std::log(alpha1);
  params[1] = std::log(alpha12);
  params[2] = std::log(alpha2);
  params[3] = rho12;

  return params;
}

unsigned int BaseLogLikelihood::GetNumberOfParameters()
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

double BaseLogLikelihood::GetLogDeterminant()
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
      double sqDist = m_DistanceMatrix(i, j) * m_DistanceMatrix(i, j);
      unsigned int workLabel = m_PointLabels[i] + m_PointLabels[j];

      if (workLabel == 2)
        resVal = this->EvaluateLFunction(sqDist, m_FirstIntensity, m_FirstAmplitude, m_FirstAlpha);
      else if (workLabel == 3)
        resVal = this->EvaluateLFunction(sqDist, m_CrossIntensity, m_CrossAmplitude, m_CrossAlpha);
      else
        resVal = this->EvaluateLFunction(sqDist, m_SecondIntensity, m_SecondAmplitude, m_SecondAlpha);

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

  m_GradientLogDeterminant.set_size(this->GetNumberOfParameters());
  m_GradientLogDeterminant[0] = arma::trace(lMatrixInverse * lMatrixDeriv1);
  m_GradientLogDeterminant[1] = arma::trace(lMatrixInverse * lMatrixDeriv2);
  m_GradientLogDeterminant[2] = arma::trace(lMatrixInverse * lMatrixDeriv3);
  m_GradientLogDeterminant[3] = arma::trace(lMatrixInverse * lMatrixDeriv4);

  return resVal;
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  bool validParams = this->CheckModelParameters(x);
  if (!validParams)
    return DBL_MAX;

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_Integral << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in evaluate");
  }

  double logLik = 2.0 * m_DomainVolume;
  logLik += m_DomainVolume * m_Integral;
  logLik += m_LogDeterminant;

  return -2.0 * logLik;
}

void BaseLogLikelihood::Gradient(const arma::mat& x, arma::mat &g)
{
  g.set_size(this->GetNumberOfParameters(), 1);

  bool validParams = this->CheckModelParameters(x);
  if (!validParams)
  {
    g.fill(0.0);
    return;
  }

  this->SetModelParameters(x);

  if (m_Modified)
  {
    m_Integral = this->GetIntegral();
    m_LogDeterminant = this->GetLogDeterminant();
  }

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_Integral << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in gradient");
  }

  for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
    g[i] = m_DomainVolume * m_GradientIntegral[i] + m_GradientLogDeterminant[i];

  g *= -2.0;
}

double BaseLogLikelihood::EvaluateWithGradient(const arma::mat& x, arma::mat& g)
{
  this->SetModelParameters(x);
  g.set_size(this->GetNumberOfParameters(), 1);

  bool validParams = this->CheckModelParameters(x);
  if (!validParams)
  {
    g.fill(0.0);
    return DBL_MAX;
  }

  m_Integral = this->GetIntegral();
  m_LogDeterminant = this->GetLogDeterminant();

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_Integral << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in evaluate with gradient");
  }

  double logLik = 2.0 * m_DomainVolume;
  logLik += m_DomainVolume * m_Integral;
  logLik += m_LogDeterminant;

  for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
    g[i] = m_DomainVolume * m_GradientIntegral[i] + m_GradientLogDeterminant[i];

  g *= -2.0;

  return -2.0 * logLik;
}

double BaseLogLikelihood::EvaluateConstraint(const size_t i, const arma::mat& x)
{
  this->CheckModelParameters(x);
  this->SetModelParameters(x);
  return m_ConstraintVector[i];
}

void BaseLogLikelihood::GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g)
{
  g.set_size(this->GetNumberOfParameters(), 1);
  g.fill(0.0);
}

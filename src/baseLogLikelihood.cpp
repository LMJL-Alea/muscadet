#include "baseLogLikelihood.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

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

arma::mat BaseLogLikelihood::GetInitialPoint()
{
  arma::mat params(this->GetNumberOfParameters(), 1);

  // Retrieve Poisson estimates of rho_1 and rho_2
  m_FirstIntensity = 0.0;
  m_SecondIntensity = 0.0;

  for (unsigned int i = 0;i < m_SampleSize;++i)
  {
    if (m_PointLabels[i] == 1)
      ++m_FirstIntensity;

    if (m_PointLabels[i] == 2)
      ++m_SecondIntensity;
  }

  m_FirstIntensity /= m_DomainVolume;
  m_SecondIntensity /= m_DomainVolume;

  bool paramsOk = false;

  while (!paramsOk)
  {
    double sample = arma::randu();

    if (m_EstimateFirstAmplitude)
      m_FirstAmplitude = sample;

    if (m_EstimateSecondAmplitude)
      m_SecondAmplitude = sample;

    unsigned int pos = 0;

    if (m_EstimateFirstAlpha)
    {
      m_FirstAlpha = this->RetrieveAlphaFromParameters(m_FirstAmplitude, m_FirstIntensity, m_DomainDimension);
      params[pos] = m_FirstAlpha;
      ++pos;
    }

    if (m_EstimateSecondAlpha)
    {
      m_SecondAlpha = this->RetrieveAlphaFromParameters(m_SecondAmplitude, m_SecondIntensity, m_DomainDimension);
      params[pos] = m_SecondAlpha;
      ++pos;
    }

    if (m_EstimateCrossAlpha)
    {
      m_CrossAlpha = std::max(m_FirstAlpha, m_SecondAlpha) + 0.01;
      params[pos] = m_CrossAlpha;
      ++pos;
    }

    if (m_EstimateFirstAmplitude)
    {
      params[pos] = m_FirstAmplitude;
      ++pos;
    }

    if (m_EstimateSecondAmplitude)
    {
      params[pos] = m_SecondAmplitude;
      ++pos;
    }

    if (m_EstimateCrossAmplitude)
    {
      m_CrossAmplitude = 0.0;
      params[pos] = m_CrossAmplitude;
    }

    paramsOk = this->CheckModelParameters();
  }

  return params;
}

unsigned int BaseLogLikelihood::GetNumberOfParameters()
{
  unsigned int numParams = 0;

  if (m_EstimateFirstAmplitude)
    ++numParams;

  if (m_EstimateSecondAmplitude)
    ++numParams;

  if (m_EstimateCrossAmplitude)
    ++numParams;

  if (m_EstimateFirstAlpha)
    ++numParams;

  if (m_EstimateSecondAlpha)
    ++numParams;

  if (m_EstimateCrossAlpha)
    ++numParams;

  return numParams;
}

double BaseLogLikelihood::GetIntegral()
{
  typedef boost::math::quadrature::gauss_kronrod<double, 15> QuadratureType;
  const double lBound = 0.0;
  const double uBound = std::numeric_limits<double>::infinity();

  BaseIntegrand integrand;
  integrand.SetKFunction(this->GetKFunction());
  integrand.SetFirstAlpha(m_FirstAlpha);
  integrand.SetSecondAlpha(m_SecondAlpha);
  integrand.SetCrossAlpha(m_CrossAlpha);
  integrand.SetFirstAmplitude(m_FirstAmplitude);
  integrand.SetSecondAmplitude(m_SecondAmplitude);
  integrand.SetCrossAmplitude(m_CrossAmplitude);
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
      double tmpVal = this->EvaluateL12Function(sqDist, m_FirstAmplitude, m_SecondAmplitude, m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);

      if (workLabel == 2)
        resVal = this->EvaluateLFunction(sqDist, m_FirstAmplitude, m_CrossAmplitude, m_FirstAlpha, m_CrossAlpha, tmpVal, m_DomainDimension);
      else if (workLabel == 3)
        resVal = tmpVal;
      else
        resVal = this->EvaluateLFunction(sqDist, m_SecondAmplitude, m_CrossAmplitude, m_SecondAlpha, m_CrossAlpha, tmpVal, m_DomainDimension);

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

  bool validParams = this->CheckModelParameters();
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
  this->SetModelParameters(x);
  g.set_size(this->GetNumberOfParameters(), 1);

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

  bool validParams = this->CheckModelParameters();
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
  this->SetModelParameters(x);
  this->CheckModelParameters();
  return m_ConstraintVector[i];
}

void BaseLogLikelihood::GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g)
{
  g.set_size(this->GetNumberOfParameters(), 1);
  g.fill(0.0);
}

void BaseLogLikelihood::SetFirstAlpha(const double x)
{
  m_FirstAlpha = x;
  m_FirstIntensity = this->RetrieveIntensityFromParameters(m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
  m_EstimateFirstAlpha = false;
}

void BaseLogLikelihood::SetSecondAlpha(const double x)
{
  m_SecondAlpha = x;
  m_SecondIntensity = this->RetrieveIntensityFromParameters(m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
  m_EstimateSecondAlpha = false;
}

void BaseLogLikelihood::SetCrossAlpha(const double x)
{
  m_CrossAlpha = x;
  m_CrossIntensity = this->RetrieveIntensityFromParameters(m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);
  m_EstimateCrossAlpha = false;
}

void BaseLogLikelihood::SetFirstAmplitude(const double x)
{
  m_FirstAmplitude = x;
  m_FirstIntensity = this->RetrieveIntensityFromParameters(m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
  m_EstimateFirstAmplitude = false;
}

void BaseLogLikelihood::SetSecondAmplitude(const double x)
{
  m_SecondAmplitude = x;
  m_SecondIntensity = this->RetrieveIntensityFromParameters(m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
  m_EstimateSecondAmplitude = false;
}

void BaseLogLikelihood::SetCrossAmplitude(const double x)
{
  m_CrossAmplitude = x;
  m_CrossIntensity = this->RetrieveIntensityFromParameters(m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);
  m_EstimateCrossAmplitude = false;
}

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  m_Modified = false;

  unsigned int pos = 0;

  if (m_EstimateFirstAlpha)
  {
    double workScalar = params[pos];

    if (m_FirstAlpha != workScalar)
    {
      m_FirstAlpha = workScalar;
      m_FirstIntensity = this->RetrieveIntensityFromParameters(m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateSecondAlpha)
  {
    double workScalar = params[pos];

    if (m_SecondAlpha != workScalar)
    {
      m_SecondAlpha = workScalar;
      m_SecondIntensity = this->RetrieveIntensityFromParameters(m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateCrossAlpha)
  {
    double workScalar = params[pos];

    if (m_CrossAlpha != workScalar)
    {
      m_CrossAlpha = workScalar;
      m_CrossIntensity = this->RetrieveIntensityFromParameters(m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateFirstAmplitude)
  {
    double workScalar = params[pos];

    if (m_FirstAmplitude != workScalar)
    {
      m_FirstAmplitude = workScalar;
      m_FirstIntensity = this->RetrieveIntensityFromParameters(m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateSecondAmplitude)
  {
    double workScalar = params[pos];

    if (m_SecondAmplitude != workScalar)
    {
      m_SecondAmplitude = workScalar;
      m_SecondIntensity = this->RetrieveIntensityFromParameters(m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  if (m_EstimateCrossAmplitude)
  {
    double workScalar = params[pos];

    if (m_CrossAmplitude != workScalar)
    {
      m_CrossAmplitude = workScalar;
      m_CrossIntensity = this->RetrieveIntensityFromParameters(m_CrossAmplitude, m_CrossAlpha, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;
  }

  Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossIntensity << " " << m_FirstAmplitude << " " << m_SecondAmplitude << " " << m_CrossAmplitude << std::endl;
}

bool BaseLogLikelihood::CheckModelParameters()
{
  if (m_FirstAmplitude > 1.0 - m_Epsilon || m_FirstAmplitude < m_Epsilon)
    return false;

  if (m_SecondAmplitude > 1.0 - m_Epsilon || m_SecondAmplitude < m_Epsilon)
    return false;

  if (m_FirstAlpha < m_Epsilon)
    return false;

  if (m_SecondAlpha < m_Epsilon)
    return false;

  if (this->EvaluateAlphaConstraint(m_FirstAlpha, m_SecondAlpha, m_CrossAlpha))
    return false;

  if (m_CrossAmplitude * m_CrossAmplitude > std::min(m_FirstAmplitude * m_SecondAmplitude, (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon))
    return false;

  return true;
}

double BaseLogLikelihood::GetBesselJRatio(const double sqDist, const double alpha, const unsigned int dimension)
{
  double order = (double)dimension / 2.0;
  double tmpVal = std::sqrt(2.0 * (double)dimension * sqDist) / alpha;

  if (tmpVal < std::sqrt(std::numeric_limits<double>::epsilon()))
    return 1.0 / boost::math::tgamma(1.0 + order);

  return boost::math::cyl_bessel_j(order, tmpVal) / std::pow(tmpVal / 2.0, order);
}

#include "baseLogLikelihood.h"

const double BaseLogLikelihood::m_Epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

double BaseLogLikelihood::GetCrossAmplitudeNormalizationFactor()
{
  double normalizationFactor = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon;
  normalizationFactor = std::min(normalizationFactor, m_FirstAmplitude * m_SecondAmplitude);
  normalizationFactor = std::max(normalizationFactor, 0.0);
  normalizationFactor = std::sqrt(normalizationFactor);
  return normalizationFactor;
}

void BaseLogLikelihood::SetFirstAlpha(const double val)
{
  m_FirstAlpha = val;
}

void BaseLogLikelihood::SetSecondAlpha(const double val)
{
  m_SecondAlpha = val;
}

void BaseLogLikelihood::SetCrossAlpha(const double val)
{
  m_CrossAlpha = val;
}

void BaseLogLikelihood::SetCorrelation(const double val)
{
  m_Correlation = val;
}

void BaseLogLikelihood::SetIntegerGrid()
{
  unsigned int kBuilderSize = 2 * (unsigned int)m_TruncationIndex + 1;
  unsigned int numberOfKVectors = std::pow(kBuilderSize, m_DomainDimension);
  m_IntegerGrid.resize(numberOfKVectors);

  if (m_UseVerbose)
    Rcpp::Rcout << "* Number of k-vectors: " << numberOfKVectors << std::endl;

  arma::ivec kBuilder = arma::linspace<arma::ivec>(-m_TruncationIndex, m_TruncationIndex, kBuilderSize);
  std::vector<int> kVector(m_DomainDimension);

  for (unsigned int i = 0;i < numberOfKVectors;++i)
  {
    auto temp = i;

    int squaredDistanceValue = 0;

    for (unsigned int j = 0;j < m_DomainDimension;++j)
    {
      auto index = temp % kBuilderSize;
      temp /= kBuilderSize;
      kVector[j] = kBuilder[index] * m_DeltaDiagonal[j];
      squaredDistanceValue += kVector[j] * kVector[j];
    }

    m_IntegerGrid[i] = std::make_pair(squaredDistanceValue, kVector);
  }

  std::sort(m_IntegerGrid.begin(), m_IntegerGrid.end());
}

void BaseLogLikelihood::SetInputData(
    const arma::mat &points,
    const arma::vec &lb,
    const arma::vec &ub,
    const arma::uvec &labels)
{
  m_DataPoints = points;
  m_DomainDimension = points.n_cols;
  m_NumberOfPoints = points.n_rows;

  m_DomainVolume = 1.0;
  m_DeltaDiagonal.set_size(m_DomainDimension);
  for (unsigned int i = 0;i < m_DomainDimension;++i)
  {
    double segmentLength = ub[i] - lb[i];
    m_DomainVolume *= segmentLength;
    m_DeltaDiagonal[i] = 1.0 / segmentLength;
  }

  m_PointLabels = labels;

  arma::uvec uniqueLabels = arma::unique(m_PointLabels);
  m_NumberOfMarks = uniqueLabels.n_elem;

  if (m_NumberOfMarks > 2)
    Rcpp::stop("The current version only handles univariate and bivariate DPPs.");

  // Poisson estimates for intensities
  m_FirstIntensity = 0.0;
  m_SecondIntensity = 0.0;

  for (unsigned int i = 0;i < m_NumberOfPoints;++i)
  {
    if (m_PointLabels[i] == 1)
      m_FirstIntensity += 1.0;
    if (m_PointLabels[i] == 2)
      m_SecondIntensity += 1.0;
  }

  m_FirstIntensity /= m_DomainVolume;
  m_SecondIntensity /= m_DomainVolume;

  // Generate K-Space grid
  this->SetIntegerGrid();

  Rcpp::Rcout << "* Number of points:             " << m_NumberOfPoints << std::endl;
  Rcpp::Rcout << "* Domain dimension:             " << m_DomainDimension << std::endl;
  Rcpp::Rcout << "* Domain volume:                " << m_DomainVolume << std::endl;
  Rcpp::Rcout << "* Number of marks:              " << m_NumberOfMarks << std::endl;
  Rcpp::Rcout << "* Intensity of the first mark:  " << m_FirstIntensity << std::endl;
  Rcpp::Rcout << "* Intensity of the second mark: " << m_SecondIntensity << std::endl;
  Rcpp::Rcout << "* Truncation index:             " << m_TruncationIndex << std::endl;
}

unsigned int BaseLogLikelihood::GetNumberOfParameters()
{
  return (m_NumberOfMarks == 1) ? 1 : 4;
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  if (m_FirstAmplitude < m_Epsilon || m_FirstAmplitude > 1.0)
    return(DBL_MAX);

  arma::vec eigenValues(m_NumberOfMarks);
  arma::mat eigenVectors(m_NumberOfMarks, m_NumberOfMarks);
  arma::mat lMatrix(m_NumberOfMarks, m_NumberOfMarks);
  m_Integral = 0.0;
  arma::mat lDataMatrix(m_NumberOfPoints, m_NumberOfPoints);
  lDataMatrix.fill(0.0);

  double traceValue = 0.0;
  unsigned counter = 0;

  for (unsigned int i = 0;i < m_IntegerGrid.size();++i)
  {
    ++counter;

    // First get current k-vector (properly scaled so in fact Delta k)
    double kSquaredNorm = m_IntegerGrid[i].first;

    // Now compute eigenvalues and vector of K0^hat(Delta k)
    double k11Value = this->GetK11Value(kSquaredNorm);
    double k12Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK12Value(kSquaredNorm);
    double k22Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK22Value(kSquaredNorm);

    traceValue += k11Value + k22Value;

    if (traceValue > 0.99 * (m_FirstIntensity + m_SecondIntensity))
      break;

    double dValue = std::sqrt((k11Value - k22Value) * (k11Value - k22Value) + 4 * k12Value * k12Value);
    eigenValues[0] = (k11Value + k22Value + dValue) / 2.0;

    if (!arma::is_finite(eigenValues[0]))
    {
      Rcpp::Rcout << k11Value << " " << k12Value << " " << k22Value << " " << dValue << " " << kSquaredNorm << std::endl;
      Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossAmplitude / M_PI / m_CrossAlpha / m_CrossAlpha / std::sqrt(m_FirstIntensity * m_SecondIntensity) << std::endl;
      Rcpp::stop("Aborting...");
    }

    if (1.0 - eigenValues[0] < 0)
      Rcpp::Rcout << "AU LOUP : " << 1.0 - eigenValues[0] << std::endl;

    if (m_NumberOfMarks == 2)
      eigenValues[1] = (k11Value + k22Value - dValue) / 2.0;

    if (m_UseVerbose)
      Rcpp::Rcout << eigenValues << std::endl;

    // Now focus on the log determinant part

    eigenVectors(0, 0) = 1.0;

    if (m_NumberOfMarks == 2)
    {
      double mValue = (-k11Value + k22Value + dValue) / (2.0 * k12Value);
      double mValueDeriv = std::sqrt(1.0 + mValue * mValue);
      eigenVectors(0, 0) = 1.0 / mValueDeriv;
      eigenVectors(1, 0) = mValue / mValueDeriv;
      eigenVectors(0, 1) = -mValue / mValueDeriv;
      eigenVectors(1, 1) = 1.0 / mValueDeriv;
    }

    if (m_UseVerbose)
      Rcpp::Rcout << eigenVectors << std::endl;

    // Compute first summation (which does not depend on either eigenvectors or data)
    // and matrix involved in log det that does not depend on data
    lMatrix.fill(0.0);
    for (unsigned int k = 0;k < m_NumberOfMarks;++k)
    {
      m_Integral += std::log(1.0 - eigenValues[k]);
      lMatrix = lMatrix + eigenValues[k] / (1.0 - eigenValues[k]) * (eigenVectors.col(k) * eigenVectors.col(k).t());
    }

    for (unsigned int l = 0;l < m_NumberOfPoints;++l)
    {
      for (unsigned int m = l;m < m_NumberOfPoints;++m)
      {
        double scalarProduct = 0.0;
        for (unsigned int j = 0;j < m_DomainDimension;++j)
          scalarProduct += (m_IntegerGrid[i].second)[j] * (m_DataPoints(l, j) - m_DataPoints(m, j));

        double workValue = std::cos(2.0 * M_PI * scalarProduct);
        workValue *= lMatrix(m_PointLabels[l] - 1, m_PointLabels[m] - 1);
        lDataMatrix(l, m) += workValue;
        if (l != m)
          lDataMatrix(m, l) += workValue;
      }
    }
  }

  if (m_UseVerbose)
    Rcpp::Rcout << "* Number of k-vector used:          " << counter << std::endl;

  double sign;
  arma::log_det(m_LogDeterminant, sign, lDataMatrix);

  if (m_UseVerbose)
    Rcpp::Rcout << "* Log-determinant: " << m_LogDeterminant << std::endl;

  if (!std::isfinite(m_Integral) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_Integral << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in evaluate");
  }

  double logLik = m_NumberOfMarks * m_DomainVolume;
  logLik += m_Integral;
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

  m_Integral = 0;
  m_LogDeterminant = 0;

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

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  // if (this->GetNumberOfParameters() == 1)
  // {
  //   m_FirstAlpha = params[0];
  //
  // }
  // else if (this->GetNumberOfParameters() == 1)
  // {
  //
  // }
  // else if (this->GetNumberOfParameters() == 1)
  // {
  //
  // }
  // else if (this->GetNumberOfParameters() == 1)
  // {
  //
  // }
  // else
  //   Rcpp::stop("The input parameter vector should be of length 1, 2, 4 or 6.");

  unsigned int pos = 0;
  double workScalar = 0.0;

  // Set k1
  workScalar = params[pos];

  if (m_FirstAmplitude != workScalar)
  {
    m_FirstAmplitude = workScalar;
    m_FirstAlpha = this->RetrieveAlphaFromParameters(m_FirstAmplitude, m_FirstIntensity, m_DomainDimension);
  }

  ++pos;

  if (m_NumberOfMarks == 2)
  {
    // Set k2
    workScalar = params[pos];

    if (m_SecondAmplitude != workScalar)
    {
      m_SecondAmplitude = workScalar;
      m_SecondAlpha = this->RetrieveAlphaFromParameters(m_SecondAmplitude, m_SecondIntensity, m_DomainDimension);
    }

    ++pos;

    // Set k12star
    workScalar = params[pos];

    if (m_NormalizedCrossAmplitude != workScalar)
    {
      m_NormalizedCrossAmplitude = workScalar;
    }

    ++pos;

    // Set beta12
    workScalar = params[pos];

    if (m_CrossBeta != workScalar)
    {
      m_CrossBeta = workScalar;
    }

    ++pos;
  }

  if (0) // (m_Modified)
  {
    if (m_NumberOfMarks == 2)
    {
      m_InverseCrossAlpha = m_CrossBeta / this->GetCrossAlphaLowerBound();
      m_CrossAlpha = 1.0 / m_InverseCrossAlpha;
      double upperBound = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude);
      upperBound = std::min(upperBound, m_FirstAmplitude * m_SecondAmplitude);
      upperBound = std::max(upperBound, 0.0);
      upperBound = std::sqrt(upperBound);
      m_CrossAmplitude = m_NormalizedCrossAmplitude * upperBound;
    }
  }

  // Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossAmplitude / M_PI / m_CrossAlpha / m_CrossAlpha / std::sqrt(m_FirstIntensity * m_SecondIntensity) << std::endl;
}

bool BaseLogLikelihood::CheckModelParameters()
{
  return true;
}

arma::mat BaseLogLikelihood::GetInitialPoint()
{
  unsigned int numParams = this->GetNumberOfParameters();
  arma::mat params(numParams, 1);

  if (numParams == 1)
  {
    m_FirstAmplitude = this->RetrieveAmplitudeFromParameters(m_FirstIntensity, m_FirstAlpha, m_DomainDimension);
    params(0, 0) = m_FirstAmplitude;
  }
  else if (numParams == 2)
  {
    m_FirstAmplitude = this->RetrieveAmplitudeFromParameters(m_FirstIntensity, m_FirstAlpha, m_DomainDimension);
    double alphaUpperBound = this->RetrieveAlphaFromParameters(1.0, 1.0 / m_DomainVolume, m_DomainDimension);
    params(0, 0) = m_FirstAmplitude;
    params(1, 0) = m_FirstAlpha / alphaUpperBound;
  }
  else if (numParams == 4)
  {
    m_FirstAmplitude = this->RetrieveAmplitudeFromParameters(m_FirstIntensity, m_FirstAlpha, m_DomainDimension);
    m_SecondAmplitude = this->RetrieveAmplitudeFromParameters(m_SecondIntensity, m_SecondAlpha, m_DomainDimension);
    m_CrossAmplitude = this->RetrieveAmplitudeFromParameters(m_Correlation * std::sqrt(m_FirstIntensity * m_SecondIntensity), m_CrossAlpha, m_DomainDimension);

    double k12UpperBound = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon;
    k12UpperBound = std::min(k12UpperBound, m_FirstAmplitude * m_SecondAmplitude);
    k12UpperBound = std::max(k12UpperBound, 0.0);
    k12UpperBound = std::sqrt(k12UpperBound);

    params(0, 0) = m_FirstAmplitude;
    params(1, 0) = m_SecondAmplitude;
    params(2, 0) = m_CrossAmplitude / k12UpperBound;
    params(3, 0) = this->GetCrossAlphaLowerBound() / m_CrossAlpha;
  }
  else if (numParams == 6)
  {
    m_FirstAmplitude = this->RetrieveAmplitudeFromParameters(m_FirstIntensity, m_FirstAlpha, m_DomainDimension);
    m_SecondAmplitude = this->RetrieveAmplitudeFromParameters(m_SecondIntensity, m_SecondAlpha, m_DomainDimension);
    m_CrossAmplitude = this->RetrieveAmplitudeFromParameters(m_Correlation * std::sqrt(m_FirstIntensity * m_SecondIntensity), m_CrossAlpha, m_DomainDimension);

    double k12UpperBound = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon;
    k12UpperBound = std::min(k12UpperBound, m_FirstAmplitude * m_SecondAmplitude);
    k12UpperBound = std::max(k12UpperBound, 0.0);
    k12UpperBound = std::sqrt(k12UpperBound);

    double alphaUpperBound = this->RetrieveAlphaFromParameters(1.0, 1.0 / m_DomainVolume, m_DomainDimension);

    params(0, 0) = m_FirstAmplitude;
    params(1, 0) = m_SecondAmplitude;
    params(2, 0) = m_CrossAmplitude / k12UpperBound;
    params(3, 0) = this->GetCrossAlphaLowerBound() / m_CrossAlpha;
    params(4, 0) = m_FirstAlpha / alphaUpperBound;
    params(5, 0) = m_SecondAlpha / alphaUpperBound;
  }
  else
    Rcpp::stop("The input vector should be of size 1, 2, 4 or 6.");

  return params;
}

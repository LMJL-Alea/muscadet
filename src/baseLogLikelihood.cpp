#include "baseLogLikelihood.h"

const double BaseLogLikelihood::m_Epsilon = 1.0e-4;

void BaseLogLikelihood::SetInputs(
    const arma::mat &points,
    const arma::vec &lb,
    const arma::vec &ub,
    const Rcpp::Nullable<arma::uvec> &labels)
{
  m_DataPoints = points;
  m_DomainDimension = points.n_cols;
  m_SampleSize = points.n_rows;

  if (labels.isNull())
  {
    m_PointLabels.set_size(m_SampleSize);
    m_PointLabels.fill(1);
    m_NumberOfMarks = 1;
  }
  else
  {
    m_PointLabels = Rcpp::as<arma::uvec>(labels);
    arma::uvec uniqueLabels = arma::unique(m_PointLabels);
    m_NumberOfMarks = uniqueLabels.n_elem;
    if (m_NumberOfMarks > 2)
      Rcpp::stop("The current version only handles univariate and bivariate DPPs.");
  }

  m_DomainVolume = 1.0;
  m_DeltaDiagonal.set_size(m_DomainDimension);
  for (unsigned int i = 0;i < m_DomainDimension;++i)
  {
    double segmentLength = ub[i] - lb[i];
    m_DomainVolume *= segmentLength;
    m_DeltaDiagonal[i] = 1.0 / segmentLength;
  }

  Rcpp::Rcout << "* Sample size: " << m_SampleSize << std::endl;
  Rcpp::Rcout << "* Domain dimension: " << m_DomainDimension << std::endl;
  Rcpp::Rcout << "* Domain volume: " << m_DomainVolume << std::endl;
}

arma::mat BaseLogLikelihood::GetInitialPoint()
{
  arma::mat params(this->GetNumberOfParameters(), 1);
  return params;
}

unsigned int BaseLogLikelihood::GetNumberOfParameters()
{
  unsigned int numParams = (m_NumberOfMarks == 1) ? 1 : 4;

  if (m_EstimateIntensities)
    numParams += m_NumberOfMarks;

  return numParams;
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  if (m_Modified)
  {
    arma::vec kVector(m_DomainDimension);
    arma::vec eigenValues(m_NumberOfMarks);
    arma::mat eigenVectors(m_NumberOfMarks, m_NumberOfMarks);
    arma::mat lMatrix(m_NumberOfMarks, m_NumberOfMarks);
    lMatrix.fill(0.0);
    m_Integral = 0.0;
    arma::mat lDataMatrix(m_SampleSize, m_SampleSize);
    lDataMatrix.fill(0.0);

    unsigned int kBuilderSize = 2 * m_TruncationIndex + 1;
    unsigned int numberOfKVectors = std::pow(kBuilderSize, m_DomainDimension);
    Rcpp::Rcout << "* Number of k-vectors: " << numberOfKVectors << std::endl;
    arma::ivec kBuilder = arma::linspace<arma::ivec>(-m_TruncationIndex, m_TruncationIndex, kBuilderSize);

    for (unsigned int i = 0;i < numberOfKVectors;++i)
    {
      // First get current k-vector (properly scaled so in fact Delta k)
      auto temp = i;
      double kSquaredNorm = 0.0;
      for (unsigned int j = 0;j < m_DomainDimension;++j)
      {
        auto index = temp % kBuilderSize;
        temp /= kBuilderSize;
        kVector[j] = kBuilder[index] * m_DeltaDiagonal[j];
        kSquaredNorm += kVector[j] * kVector[j];
      }

      // Now compute eigenvalues and vector of K0^hat(Delta k)
      double k11Value = this->GetK11Value(kSquaredNorm);
      double k12Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK12Value(kSquaredNorm);
      double k22Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK22Value(kSquaredNorm);

      double dValue = std::sqrt((k11Value - k22Value) * (k11Value - k22Value) + 4 * k12Value * k12Value);
      eigenValues[0] = (k11Value + k22Value + dValue) / 2.0;

      if (m_NumberOfMarks == 2)
        eigenValues[1] = (k11Value + k22Value - dValue) / 2.0;

      // Now focus on the log determinant part

      double z1Value = std::sqrt(2.0 * dValue * (k22Value - k11Value + dValue));
      eigenVectors(0, 0) = (m_NumberOfMarks == 1) ? 1.0 : 2.0 * k12Value / z1Value;

      if (m_NumberOfMarks == 2)
      {
        eigenVectors(1, 0) = (k22Value - k11Value + dValue) / z1Value;
        double z2Value = std::sqrt(2.0 * dValue * (k11Value - k22Value + dValue));
        eigenVectors(0, 1) = 2.0 * k12Value / z2Value;
        eigenVectors(1, 1) = (k22Value - k11Value - dValue) / z2Value;
      }

      // Compute first summation (which does not depend on either eigenvectors or data)
      // and matrix involved in log det that does not depend on data
      for (unsigned int k = 0;k < m_NumberOfMarks;++k)
      {
        m_Integral += std::log(1.0 - eigenValues[k]);
        lMatrix = lMatrix + eigenValues[k] / (1.0 - eigenValues[k]) * (eigenVectors.col(k) * eigenVectors.col(k).t());
      }

      for (unsigned int l = 0;l < m_SampleSize;++l)
      {
        for (unsigned int m = l;m < m_SampleSize;++m)
        {
          double scalarProduct = 0.0;
          for (unsigned int k = 0;k < m_DomainDimension;++k)
            scalarProduct += kVector[k] * (m_DataPoints(l, k) - m_DataPoints(m, k));

          double workValue = std::cos(2.0 * M_PI * scalarProduct);
          workValue *= lMatrix(m_PointLabels[l] - 1, m_PointLabels[m] - 1);
          lDataMatrix(l, m) += workValue;
          if (l != m)
            lDataMatrix(m, l) += workValue;
        }
      }
    }

    // std::complex<double> logDeterminant = arma::log_det(lDataMatrix);
    // Rcpp::Rcout << "* Complex log-determinant: " << logDeterminant << std::endl;
    // m_LogDeterminant = std::real(logDeterminant);

    double sign;
    arma::log_det(m_LogDeterminant, sign, lDataMatrix);
    Rcpp::Rcout << "* Complex log-determinant: " << m_LogDeterminant << std::endl;
  }

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

  if (m_Modified)
  {
    m_Integral = 0;
    m_LogDeterminant = 0;
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

void BaseLogLikelihood::SetIntensities(const double rho1, const double rho2)
{
  m_FirstIntensity = rho1;

  if (m_NumberOfMarks == 1)
  {
    m_EstimateIntensities = false;
    return;
  }

  if (rho2 == NA_REAL)
    Rcpp::stop("For bivariate models, you need to supply a value for the second marginal intensity as well.");

  m_SecondIntensity = rho2;
  m_EstimateIntensities = false;
}

void BaseLogLikelihood::SetModelParameters(const arma::mat &params)
{
  if (m_Modified)
    return;

  unsigned int pos = 0;
  double workScalar = 0.0;

  // Set k1
  workScalar = params[pos];

  if (m_FirstAmplitude != workScalar)
  {
    m_FirstAmplitude = workScalar;
    if (!m_EstimateIntensities)
      m_FirstAlpha = this->RetrieveAlphaFromParameters(m_FirstAmplitude, m_FirstIntensity, m_DomainDimension);
    m_Modified = true;
  }

  ++pos;

  if (m_NumberOfMarks == 2)
  {
    // Set k2
    workScalar = params[pos];

    if (m_SecondAmplitude != workScalar)
    {
      m_SecondAmplitude = workScalar;
      if (!m_EstimateIntensities)
        m_SecondAlpha = this->RetrieveAlphaFromParameters(m_SecondAmplitude, m_SecondIntensity, m_DomainDimension);
      m_Modified = true;
    }

    ++pos;

    // Set k12star
    workScalar = params[pos];

    if (m_NormalizedCrossAmplitude != workScalar)
    {
      m_NormalizedCrossAmplitude = workScalar;
      m_Modified = true;
    }

    ++pos;

    // Set beta12
    workScalar = params[pos];

    if (m_CrossBeta != workScalar)
    {
      m_CrossBeta = workScalar;
      m_Modified = true;
    }

    ++pos;
  }

  // Set alpha_i_star
  if (m_EstimateIntensities)
  {
    double upperBound = this->RetrieveAlphaFromParameters(1.0, 1.0 / m_DomainVolume, m_DomainDimension);

    workScalar = params[pos];

    if (m_NormalizedFirstAlpha != workScalar)
    {
      m_NormalizedFirstAlpha = workScalar;
      m_FirstAlpha = m_NormalizedFirstAlpha * upperBound;
      m_Modified = true;
    }

    ++pos;

    if (m_NumberOfMarks == 2)
    {
      workScalar = params[pos];

      if (m_NormalizedSecondAlpha != workScalar)
      {
        m_NormalizedSecondAlpha = workScalar;
        m_SecondAlpha = m_NormalizedSecondAlpha * upperBound;
        m_Modified = true;
      }
    }
  }

  if (m_Modified)
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

    if (m_EstimateIntensities)
    {
      m_FirstIntensity = this->RetrieveIntensityFromParameters(m_FirstAmplitude, m_FirstAlpha, m_DomainDimension);
      if (m_NumberOfMarks == 2)
        m_SecondIntensity = this->RetrieveIntensityFromParameters(m_SecondAmplitude, m_SecondAlpha, m_DomainDimension);
    }
  }

  // Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossAmplitude / M_PI / m_CrossAlpha / m_CrossAlpha / std::sqrt(m_FirstIntensity * m_SecondIntensity) << std::endl;
}

bool BaseLogLikelihood::CheckModelParameters()
{
  return true;
}

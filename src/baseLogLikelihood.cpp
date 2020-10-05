#include "baseLogLikelihood.h"

const double BaseLogLikelihood::m_Epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

void BaseLogLikelihood::SetFirstAlpha(const double val)
{
  m_FirstAlpha = val;
  m_FirstAmplitude = this->RetrieveAmplitudeFromParameters(m_FirstIntensity, val, m_DomainDimension);
}

void BaseLogLikelihood::SetSecondAlpha(const double val)
{
  m_SecondAlpha = val;
  m_SecondAmplitude = this->RetrieveAmplitudeFromParameters(m_SecondIntensity, val, m_DomainDimension);
}

double BaseLogLikelihood::GetCrossAmplitudeNormalizationFactor()
{
  double normalizationFactor = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon;
  normalizationFactor = std::min(normalizationFactor, m_FirstAmplitude * m_SecondAmplitude);
  normalizationFactor = std::max(normalizationFactor, 0.0);
  normalizationFactor = std::sqrt(normalizationFactor);
  return normalizationFactor;
}

void BaseLogLikelihood::SetInputData(
    const arma::mat &points,
    const arma::vec &lb,
    const arma::vec &ub,
    const arma::uvec &labels,
    const Rcpp::DataFrame &ndGrid,
    const unsigned int N)
{
  // m_DataPoints = points;
  m_TruncationIndex = N;
  m_DomainDimension = points.n_cols;
  m_NumberOfPoints = points.n_rows;
  m_DataLMatrix.set_size(m_NumberOfPoints, m_NumberOfPoints);

  m_CosineMatrix.set_size(m_NumberOfPoints * (m_NumberOfPoints - 1) / 2, m_DomainDimension);
  unsigned int pos = 0;
  for (unsigned int i = 0;i < m_NumberOfPoints;++i)
  {
    for (unsigned int j = i + 1;j < m_NumberOfPoints;++j)
    {
      m_CosineMatrix.row(pos) = arma::cos(2.0 * M_PI * (points.row(i) - points.row(j)));
      ++pos;
    }
  }

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

  m_InternalLMatrix.set_size(m_NumberOfMarks, m_NumberOfMarks);
  m_WorkingEigenValues.set_size(m_NumberOfMarks);
  m_WorkingEigenVectors.set_size(m_NumberOfMarks, m_NumberOfMarks);
  m_LMatrixSum.set_size(m_NumberOfMarks, m_NumberOfMarks);

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
  m_MaximalNumberOfKVectors = ndGrid.nrows();
  m_KGrid = static_cast<Rcpp::IntegerMatrix>(Rcpp::no_init(m_MaximalNumberOfKVectors, m_DomainDimension));
  for (unsigned int i = 0;i < m_DomainDimension;++i)
    m_KGrid.column(i) = Rcpp::as<Rcpp::IntegerVector>(ndGrid[i]);
  m_KSquaredNorms = ndGrid["sq_norm"];
  m_KWeights = ndGrid["weight"];

  m_ListOfInternalLMatrices.set_size(m_NumberOfMarks, m_NumberOfMarks, m_MaximalNumberOfKVectors);
  m_CosineValues.set_size(m_TruncationIndex + 1, m_DomainDimension);
  m_CosineValues.row(0).fill(1.0);

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

void BaseLogLikelihood::ComputeAll()
{
  m_LogSpectrum = 0.0;
  double traceValue = 0.0;
  double kSquaredNorm = 0.0;
  double prevKSquaredNorm = 0.0;
  m_DataLMatrix.fill(0.0);
  double sign = 0.0;

  for (unsigned int i = 0;i < m_MaximalNumberOfKVectors;++i)
  {
    kSquaredNorm = m_KSquaredNorms[i];
    double weightValue = m_KWeights[i];

    // Now compute eigenvalues and vector of K0^hat(Delta k)
    double k11Value = this->GetK11Value(kSquaredNorm);
    double k22Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK22Value(kSquaredNorm);

    traceValue += (k11Value + k22Value) * weightValue;

    // Rcpp::Rcout << i << " " << traceValue << " " << kSquaredNorm << std::endl;

    if (traceValue > m_RelativeTolerance * (m_FirstIntensity + m_SecondIntensity) && kSquaredNorm != prevKSquaredNorm)
      break;

    double k12Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK12Value(kSquaredNorm);

    double dValue = std::sqrt((k11Value - k22Value) * (k11Value - k22Value) + 4.0 * k12Value * k12Value);
    m_WorkingEigenValues[0] = (k11Value + k22Value + dValue) / 2.0;

    if (m_NumberOfMarks == 2)
      m_WorkingEigenValues[1] = (k11Value + k22Value - dValue) / 2.0;

    // Now focus on the log determinant part

    m_WorkingEigenVectors(0, 0) = 1.0;

    if (m_NumberOfMarks == 2)
    {
      double mValue = (-k11Value + k22Value + dValue) / (2.0 * k12Value);
      double mValueDeriv = std::sqrt(1.0 + mValue * mValue);
      m_WorkingEigenVectors(0, 0) = 1.0 / mValueDeriv;
      m_WorkingEigenVectors(1, 0) = mValue / mValueDeriv;
      m_WorkingEigenVectors(0, 1) = -mValue / mValueDeriv;
      m_WorkingEigenVectors(1, 1) = 1.0 / mValueDeriv;
    }

    // Compute first summation (which does not depend on either eigenvectors or data)
    // and matrix involved in log det that does not depend on data
    m_InternalLMatrix.fill(0.0);
    for (unsigned int j = 0;j < m_NumberOfMarks;++j)
    {
      m_LogSpectrum += std::log(1.0 - m_WorkingEigenValues[j]) * weightValue;
      m_InternalLMatrix = m_InternalLMatrix + m_WorkingEigenValues[j] / (1.0 - m_WorkingEigenValues[j]) * (m_WorkingEigenVectors.col(j) * m_WorkingEigenVectors.col(j).t());
    }
    m_InternalLMatrix = m_InternalLMatrix * weightValue;

    unsigned int pos = 0;
    unsigned int nValue = 0;
    for (unsigned int j = 0;j < m_DomainDimension;++j)
      nValue = std::max(nValue, (unsigned int)m_KGrid(i, j));

    for (unsigned int j = 0;j < m_NumberOfPoints;++j)
    {
      m_DataLMatrix(j, j) += m_InternalLMatrix(m_PointLabels[j] - 1, m_PointLabels[j] - 1);

      for (unsigned int k = j + 1;k < m_NumberOfPoints;++k)
      {
        m_CosineValues.row(1) = m_CosineMatrix.row(pos);
        for (unsigned int l = 2;l < nValue + 1;++l)
          m_CosineValues.row(l) = 2.0 * m_CosineValues.row(1) % m_CosineValues.row(l - 1) - m_CosineValues.row(l - 2);

        double workValue = m_InternalLMatrix(m_PointLabels[j] - 1, m_PointLabels[k] - 1);
        for (unsigned int l = 0;l < m_DomainDimension;++l)
          workValue *= m_CosineValues(m_KGrid(i, l), l);

        m_DataLMatrix(j, k) += workValue;
        m_DataLMatrix(k, j) += workValue;

        ++pos;
      }
    }

    prevKSquaredNorm = kSquaredNorm;
  }

  m_DataLMatrix.clean(arma::datum::eps);
  arma::log_det(m_LogDeterminant, sign, m_DataLMatrix);
}

void BaseLogLikelihood::ComputeLogSpectrum()
{
  m_LogSpectrum = 0.0;
  double traceValue = 0.0;
  m_NumberOfKVectors = 0;
  m_ActualTruncationIndex = 0;
  m_LMatrixSum.fill(0.0);
  double prevKSquaredNorm = 0.0;

  for (unsigned int i = 0;i < m_MaximalNumberOfKVectors;++i)
  {
    double kSquaredNorm = m_KSquaredNorms[i];
    double weightValue = m_KWeights[i];

    // Now compute eigenvalues and vector of K0^hat(Delta k)
    double k11Value = this->GetK11Value(kSquaredNorm);
    double k22Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK22Value(kSquaredNorm);

    traceValue += (k11Value + k22Value) * weightValue;

    // Rcpp::Rcout << i << " " << traceValue << " " << kSquaredNorm << std::endl;

    if (traceValue > m_RelativeTolerance * (m_FirstIntensity + m_SecondIntensity) && kSquaredNorm != prevKSquaredNorm)
      break;

    double k12Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK12Value(kSquaredNorm);

    for (unsigned int j = 0;j < m_DomainDimension;++j)
      m_ActualTruncationIndex = std::max(m_ActualTruncationIndex, (unsigned int)m_KGrid(i, j));

    double dValue = std::sqrt((k11Value - k22Value) * (k11Value - k22Value) + 4.0 * k12Value * k12Value);
    m_WorkingEigenValues[0] = (k11Value + k22Value + dValue) / 2.0;

    if (m_NumberOfMarks == 2)
      m_WorkingEigenValues[1] = (k11Value + k22Value - dValue) / 2.0;

    // Now focus on the log determinant part

    m_WorkingEigenVectors(0, 0) = 1.0;

    if (m_NumberOfMarks == 2)
    {
      double mValue = (-k11Value + k22Value + dValue) / (2.0 * k12Value);
      double mValueDeriv = std::sqrt(1.0 + mValue * mValue);
      m_WorkingEigenVectors(0, 0) = 1.0 / mValueDeriv;
      m_WorkingEigenVectors(1, 0) = mValue / mValueDeriv;
      m_WorkingEigenVectors(0, 1) = -mValue / mValueDeriv;
      m_WorkingEigenVectors(1, 1) = 1.0 / mValueDeriv;
    }

    // Compute first summation (which does not depend on either eigenvectors or data)
    // and matrix involved in log det that does not depend on data
    m_InternalLMatrix.fill(0.0);
    for (unsigned int j = 0;j < m_NumberOfMarks;++j)
    {
      m_LogSpectrum += std::log(1.0 - m_WorkingEigenValues[j]) * weightValue;
      m_InternalLMatrix = m_InternalLMatrix + m_WorkingEigenValues[j] / (1.0 - m_WorkingEigenValues[j]) * (m_WorkingEigenVectors.col(j) * m_WorkingEigenVectors.col(j).t());
    }
    m_InternalLMatrix *= weightValue;

    m_LMatrixSum = m_LMatrixSum + m_InternalLMatrix;
    m_ListOfInternalLMatrices.slice(m_NumberOfKVectors) = m_InternalLMatrix;

    prevKSquaredNorm = kSquaredNorm;
    ++m_NumberOfKVectors;
  }

  if (m_UseVerbose)
    Rcpp::Rcout << "* Current truncation index:     " << m_ActualTruncationIndex << std::endl;
}

void BaseLogLikelihood::ComputeLogDeterminant()
{
  m_DataLMatrix.fill(0.0);
  double sign = 0.0;
  unsigned int pos = 0;

  for (unsigned int i = 0;i < m_NumberOfPoints;++i)
  {
    m_DataLMatrix(i, i) = m_LMatrixSum(m_PointLabels[i] - 1, m_PointLabels[i] - 1);

    for (unsigned int j = i + 1;j < m_NumberOfPoints;++j)
    {
      m_CosineValues.row(1) = m_CosineMatrix.row(pos);
      for (unsigned int k = 2;k < m_ActualTruncationIndex + 1;++k)
        m_CosineValues.row(k) = 2.0 * m_CosineValues.row(1) % m_CosineValues.row(k - 1) - m_CosineValues.row(k - 2);

      for (unsigned int k = 0;k < m_NumberOfKVectors;++k)
      {
        double workValue = m_ListOfInternalLMatrices(m_PointLabels[i] - 1, m_PointLabels[j] - 1, k);
        for (unsigned int l = 0;l < m_DomainDimension;++l)
          workValue *= m_CosineValues(m_KGrid(k, l), l);

        m_DataLMatrix(i, j) += workValue;
        m_DataLMatrix(j, i) += workValue;
      }

      ++pos;
    }
  }

  m_DataLMatrix.clean(arma::datum::eps);
  arma::log_det(m_LogDeterminant, sign, m_DataLMatrix);
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  if (m_FirstAmplitude < m_Epsilon || m_FirstAmplitude > 1.0)
    return(DBL_MAX);

  // this->ComputeAll();

  this->ComputeLogSpectrum();
  this->ComputeLogDeterminant();

  if (m_UseVerbose)
  {
    Rcpp::Rcout << "* Log-spectrum: " << m_LogSpectrum << std::endl;
    Rcpp::Rcout << "* Log-determinant: " << m_LogDeterminant << std::endl;
  }

  if (!std::isfinite(m_LogSpectrum) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_LogSpectrum << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in evaluate");
  }

  double logLik = m_NumberOfMarks * m_DomainVolume;
  logLik += m_LogSpectrum;
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

  if (!std::isfinite(m_LogSpectrum) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_LogSpectrum << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
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

  m_LogSpectrum = 0;
  m_LogDeterminant = 0;

  if (!std::isfinite(m_LogSpectrum) || !std::isfinite(m_LogDeterminant))
  {
    Rcpp::Rcout << m_LogSpectrum << " " << m_LogDeterminant << " " << x.as_row() << std::endl;
    Rcpp::stop("Non finite stuff in evaluate with gradient");
  }

  double logLik = 2.0 * m_DomainVolume;
  logLik += m_DomainVolume * m_LogSpectrum;
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

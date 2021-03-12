#include "baseLogLikelihood.h"

#ifdef _OPENMP
#include <omp.h>
#endif

const double BaseLogLikelihood::m_Epsilon = std::sqrt(arma::datum::eps);

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

void BaseLogLikelihood::SetCrossParameters(const double alpha12, const double tau)
{
  m_Correlation = tau;

  if (m_Correlation < m_Epsilon)
  {
    m_Correlation = 0.0;
    m_CrossAlpha = std::sqrt(this->GetSquaredCrossAlphaLowerBound());
    m_CrossAmplitude = 0.0;
    return;
  }

  m_CrossAlpha = alpha12;
  m_CrossAmplitude = this->RetrieveAmplitudeFromParameters(
    tau * std::sqrt(m_FirstIntensity * m_SecondIntensity),
    alpha12,
    m_DomainDimension
  );
}

double BaseLogLikelihood::GetSquaredCrossAlphaUpperBound() const
{
  double weightValue = std::min(1.0, (1.0 / m_FirstAmplitude - 1.0) * (1.0/ m_SecondAmplitude - 1.0));
  double squaredCorrelationUpperBound = weightValue * std::pow(m_FirstAlpha * m_SecondAlpha / this->GetSquaredCrossAlphaLowerBound(), (double)m_DomainDimension);
  double squaredCorrelationLowerBound = 0.1 * squaredCorrelationUpperBound;
  return m_FirstAlpha * m_SecondAlpha * std::pow(weightValue / squaredCorrelationLowerBound, 1.0 / m_DomainDimension);
}

double BaseLogLikelihood::GetSquaredCrossAmplitudeUpperBound() const
{
  double normalizationFactor = (1.0 - m_FirstAmplitude) * (1.0 - m_SecondAmplitude) - m_Epsilon * m_Epsilon;
  normalizationFactor = std::min(normalizationFactor, m_FirstAmplitude * m_SecondAmplitude);
  normalizationFactor = std::max(normalizationFactor, 0.0);
  return normalizationFactor;
}

void BaseLogLikelihood::SetInputData(
    const arma::mat &points,
    const arma::vec &lb,
    const arma::vec &ub,
    const arma::uvec &labels,
    const Rcpp::DataFrame &ndGrid,
    const unsigned int N,
    const Rcpp::Nullable<arma::vec> &marginal_parameters)
{
  m_DataPoints = points;
  m_TruncationIndex = N;
  m_DomainDimension = points.n_cols;
  m_NumberOfPoints = points.n_rows;
  m_DataLMatrix.set_size(m_NumberOfPoints, m_NumberOfPoints);

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

  if (marginal_parameters.isNotNull())
  {
    // Set marginal parameters
    arma::vec marginalParams = Rcpp::as<arma::vec>(marginal_parameters);
    this->SetFirstAlpha(marginalParams(0));
    this->SetSecondAlpha(marginalParams(1));
    m_UseFixedMarginalParameters = true;
  }

  m_NumberOfParameters = (m_NumberOfMarks == 1) ? 1 : ((m_UseFixedMarginalParameters) ? 2 : 4);

  // Generate K-Space grid
  m_MaximalNumberOfKVectors = ndGrid.nrows();
  m_KGrid = static_cast<Rcpp::IntegerMatrix>(Rcpp::no_init(m_MaximalNumberOfKVectors, m_DomainDimension));
  for (unsigned int i = 0;i < m_DomainDimension;++i)
    m_KGrid.column(i) = Rcpp::as<Rcpp::IntegerVector>(ndGrid[i]);
  m_KSquaredNorms = ndGrid["sq_norm"];
  m_KWeights = ndGrid["weight"];

  m_ListOfInternalLMatrices.set_size(m_NumberOfMarks, m_NumberOfMarks, m_MaximalNumberOfKVectors);

  if (m_VerboseLevel == 1)
  {
    Rcpp::Rcout << "* Number of points:             " << m_NumberOfPoints << std::endl;
    Rcpp::Rcout << "* Domain dimension:             " << m_DomainDimension << std::endl;
    Rcpp::Rcout << "* Domain volume:                " << m_DomainVolume << std::endl;
    Rcpp::Rcout << "* Number of marks:              " << m_NumberOfMarks << std::endl;
    Rcpp::Rcout << "* Intensity of the first mark:  " << m_FirstIntensity << std::endl;
    Rcpp::Rcout << "* Intensity of the second mark: " << m_SecondIntensity << std::endl;
    Rcpp::Rcout << "* Truncation index:             " << m_TruncationIndex << std::endl;
  }
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
      double mValue = (k12Value < m_Epsilon) ? 0.0 : (-k11Value + k22Value + dValue) / (2.0 * k12Value);
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

  if (m_VerboseLevel > 1)
    Rcpp::Rcout << "* Current truncation index:     " << m_ActualTruncationIndex << std::endl;
}

void BaseLogLikelihood::ComputeLogDeterminant()
{
  m_DataLMatrix.fill(0.0);
  double sign = 0.0;
  double workInteger = 2.0 * m_NumberOfPoints + 1.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

  for (unsigned int k = 0;k < m_NumberOfPoints * (m_NumberOfPoints + 1) / 2;++k)
  {
    double i = std::floor((workInteger - std::sqrt(workInteger * workInteger - 8.0 * k)) / 2.0);
    double j = k + i * (i - 1) / 2 - (m_NumberOfPoints - 1.0) * i;
    i = (unsigned int)i;
    j = (unsigned int)j;

    if (i == j)
    {
      m_DataLMatrix(i, i) = m_LMatrixSum(m_PointLabels[i] - 1, m_PointLabels[i] - 1);
      continue;
    }

    arma::mat cosineValues(m_ActualTruncationIndex + 1, m_DomainDimension);
    cosineValues.row(0).fill(1.0);
    cosineValues.row(1) = arma::cos(2.0 * M_PI * (m_DataPoints.row(i) - m_DataPoints.row(j)));
    for (unsigned int l = 2;l < m_ActualTruncationIndex + 1;++l)
      cosineValues.row(l) = 2.0 * cosineValues.row(1) % cosineValues.row(l - 1) - cosineValues.row(l - 2);

    for (unsigned int l = 0;l < m_NumberOfKVectors;++l)
    {
      double workValue = m_ListOfInternalLMatrices(m_PointLabels[i] - 1, m_PointLabels[j] - 1, l);
      for (unsigned int m = 0;m < m_DomainDimension;++m)
        workValue *= cosineValues(m_KGrid(l, m), m);

      m_DataLMatrix(i, j) += workValue;
      m_DataLMatrix(j, i) += workValue;
    }
  }

  m_DataLMatrix.clean(arma::datum::eps);
  arma::log_det(m_LogDeterminant, sign, m_DataLMatrix);
}

double BaseLogLikelihood::GetValue(const arma::vec& x)
{
  this->SetModelParameters(x);

  if (!this->CheckModelParameters())
    return DBL_MAX;

  this->ComputeLogSpectrum();
  this->ComputeLogDeterminant();

  if (m_VerboseLevel > 1)
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

void BaseLogLikelihood::SetModelParameters(const arma::vec &params)
{
  unsigned int numParams = params.n_rows;

  if (numParams != m_NumberOfParameters)
    Rcpp::stop("The number of input parameters does not match the expected value.");

  if (numParams == 1)
    this->SetFirstAlpha(params(0));
  else if (numParams == 2)
    this->SetCrossParameters(params(0), params(1));
  else if (numParams == 4)
  {
    this->SetFirstAlpha(params(0));
    this->SetSecondAlpha(params(1));
    this->SetCrossParameters(params(2), params(3));
  }
  else
    Rcpp::stop("Unknowm number of parameters.");
}

bool BaseLogLikelihood::CheckModelParameters()
{
  if (m_FirstAmplitude < m_Epsilon || m_FirstAmplitude > 1.0 - m_Epsilon)
    return false;

  if (m_NumberOfParameters == 2)
  {
    if (m_CrossAlpha * m_CrossAlpha < this->GetSquaredCrossAlphaLowerBound())
      return false;

    if (m_CrossAmplitude * m_CrossAmplitude > this->GetSquaredCrossAmplitudeUpperBound())
      return false;
  }

  if (m_NumberOfParameters == 4)
  {
    if (m_SecondAmplitude < m_Epsilon || m_SecondAmplitude > 1.0 - m_Epsilon)
      return false;

    if (m_CrossAlpha * m_CrossAlpha < this->GetSquaredCrossAlphaLowerBound())
      return false;

    if (m_CrossAmplitude * m_CrossAmplitude > this->GetSquaredCrossAmplitudeUpperBound())
      return false;
  }

  return true;
}

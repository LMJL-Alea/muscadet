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

void BaseLogLikelihood::GenerateCombinations(const unsigned int N, const unsigned int K, std::vector<std::vector<unsigned int> > &resVector)
{
  resVector.clear();

  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's
  std::vector<unsigned int> singleCombination(K);

  // print integers and permute bitmask
  do {
    unsigned int pos = 0;
    for (unsigned int i = 0;i < N;++i) // [0..N-1] integers
    {
      if (bitmask[i])
      {
        singleCombination[pos] = i;
        ++pos;
      }
    }
    resVector.push_back(singleCombination);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

void BaseLogLikelihood::SetIntegerGrid()
{
  KVectorType kVector(m_DomainDimension);
  m_OptimizedIntegerGrid.resize(m_DomainDimension + 1);
  std::vector<std::vector<unsigned int> > combinationVector;

  // Central point
  unsigned int numberOfElements = 1;
  m_OptimizedIntegerGrid[0].resize(numberOfElements);
  std::fill(kVector.begin(), kVector.end(), 0.0);
  m_OptimizedIntegerGrid[0][0] = std::make_pair(0.0, kVector);

  // Number of points on domain axes
  this->GenerateCombinations(m_DomainDimension, 1, combinationVector);
  numberOfElements = combinationVector.size() * std::pow(m_TruncationIndex, 1);
  m_OptimizedIntegerGrid[1].resize(numberOfElements);
  for (unsigned int i = 0;i < combinationVector.size();++i)
  {
    std::fill(kVector.begin(), kVector.end(), 0.0);
    unsigned int index = combinationVector[i][0];
    for (unsigned int j = 0;j < m_TruncationIndex;++j)
    {
      kVector[index] = (j + 1.0) * m_DeltaDiagonal[index];
      m_OptimizedIntegerGrid[1][i * m_TruncationIndex + j] = std::make_pair(kVector[index] * kVector[index], kVector);
    }
  }

  // Number of points on domain planes
  this->GenerateCombinations(m_DomainDimension, 2, combinationVector);
  numberOfElements = combinationVector.size() * std::pow(m_TruncationIndex, 2);
  m_OptimizedIntegerGrid[2].resize(numberOfElements);
  unsigned int pos = 0;
  for (unsigned int i = 0;i < combinationVector.size();++i)
  {
    std::fill(kVector.begin(), kVector.end(), 0.0);
    unsigned int index1 = combinationVector[i][0];
    unsigned int index2 = combinationVector[i][1];
    for (unsigned int j = 0;j < m_TruncationIndex;++j)
    {
      kVector[index1] = (j + 1.0) * m_DeltaDiagonal[index1];
      for (unsigned int k = 0;k < m_TruncationIndex;++k)
      {
        kVector[index2] = (k + 1.0) * m_DeltaDiagonal[index2];
        m_OptimizedIntegerGrid[2][pos] = std::make_pair(kVector[index1] * kVector[index1] + kVector[index2] * kVector[index2], kVector);
        ++pos;
      }
    }
  }

  for (unsigned int i = 0;i <= m_DomainDimension;++i)
    std::sort(m_OptimizedIntegerGrid[i].begin(), m_OptimizedIntegerGrid[i].end());

  unsigned int kBuilderSize = 2 * m_TruncationIndex + 1;
  unsigned int numberOfKVectors = std::pow(kBuilderSize, m_DomainDimension);
  m_IntegerGrid.resize(numberOfKVectors);

  if (m_UseVerbose)
    Rcpp::Rcout << "* Number of k-vectors: " << numberOfKVectors << std::endl;

  arma::ivec kBuilder = arma::linspace<arma::ivec>(-m_TruncationIndex, m_TruncationIndex, kBuilderSize);

  for (unsigned int i = 0;i < numberOfKVectors;++i)
  {
    auto temp = i;

    double squaredDistanceValue = 0.0;

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

  // Rcpp::Rcout << "Integer grid: " << std::endl;
  // for (unsigned int i = 0;i < m_IntegerGrid.size();++i)
  //   Rcpp::Rcout << m_IntegerGrid[i].second[0] << " " << m_IntegerGrid[i].second[1] << std::endl;
  //
  // Rcpp::Rcout << "------------------------" << std::endl;
  // Rcpp::Rcout << "Optimized Integer grid: " << std::endl;
  // for (unsigned int i = 0;i < m_OptimizedIntegerGrid.size();++i)
  // {
  //   Rcpp::Rcout << "Dimension: " << i << std::endl;
  //   for (unsigned int j = 0;j < m_OptimizedIntegerGrid[i].size();++j)
  //     Rcpp::Rcout << m_OptimizedIntegerGrid[i][j].second[0] << " " << m_OptimizedIntegerGrid[i][j].second[1] << std::endl;
  // }
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

// double BaseLogLikelihood::GetCosN(const double x, const unsigned int n)
// {
//   if (n == 0)
//     return 1.0;
//
//   if (n == 1)
//     return std::cos((double)x);
//
//   if (n == 2)
//     return std::cos(2.0 * x);
//
//   return 2.0 * this->GetCosN(x, 1) * this->GetCosN(x, n - 1) - this->  GetCosN(x, n - 2);
// }

void BaseLogLikelihood::IncrementSummation(const double kSquaredNorm, const KVectorType &kVector, const double weight)
{
  // Now compute eigenvalues and vector of K0^hat(Delta k)
  double k11Value = this->GetK11Value(kSquaredNorm);
  double k12Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK12Value(kSquaredNorm);
  double k22Value = (m_NumberOfMarks == 1) ? 0.0 : this->GetK22Value(kSquaredNorm);

  m_TraceValue += (k11Value + k22Value) * weight;

  if (m_TraceValue > m_RelativeTolerance * (m_FirstIntensity + m_SecondIntensity))
  {
    m_ContinueLoop = false;
    return;
  }

  double dValue = std::sqrt((k11Value - k22Value) * (k11Value - k22Value) + 4 * k12Value * k12Value);
  m_WorkingEigenValues[0] = (k11Value + k22Value + dValue) / 2.0;

  if (!arma::is_finite(m_WorkingEigenValues[0]))
  {
    Rcpp::Rcout << k11Value << " " << k12Value << " " << k22Value << " " << dValue << " " << kSquaredNorm << std::endl;
    Rcpp::Rcout << m_FirstAlpha << " " << m_SecondAlpha << " " << m_CrossAlpha << " " << m_FirstIntensity << " " << m_SecondIntensity << " " << m_CrossAmplitude / M_PI / m_CrossAlpha / m_CrossAlpha / std::sqrt(m_FirstIntensity * m_SecondIntensity) << std::endl;
    Rcpp::stop("Aborting...");
  }

  if (1.0 - m_WorkingEigenValues[0] < 0)
    Rcpp::Rcout << "AU LOUP : " << 1.0 - m_WorkingEigenValues[0] << std::endl;

  if (m_NumberOfMarks == 2)
    m_WorkingEigenValues[1] = (k11Value + k22Value - dValue) / 2.0;

  if (m_UseVerbose)
    Rcpp::Rcout << m_WorkingEigenValues << std::endl;

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

  if (m_UseVerbose)
    Rcpp::Rcout << m_WorkingEigenVectors << std::endl;

  // Compute first summation (which does not depend on either eigenvectors or data)
  // and matrix involved in log det that does not depend on data
  m_InternalLMatrix.fill(0.0);
  for (unsigned int k = 0;k < m_NumberOfMarks;++k)
  {
    m_LogSpectrum += std::log(1.0 - m_WorkingEigenValues[k]) * weight;
    m_InternalLMatrix = m_InternalLMatrix + m_WorkingEigenValues[k] / (1.0 - m_WorkingEigenValues[k]) * (m_WorkingEigenVectors.col(k) * m_WorkingEigenVectors.col(k).t());
  }
  m_InternalLMatrix = m_InternalLMatrix * weight;

  double workValue = 0.0;

  for (unsigned int l = 0;l < m_NumberOfPoints;++l)
  {
    for (unsigned int m = l;m < m_NumberOfPoints;++m)
    {
      if (!m_SplitSummation)
      {
        double scalarProduct = 0.0;
        for (unsigned int k = 0;k < m_DomainDimension;++k)
          scalarProduct += kVector[k] * (m_DataPoints(l, k) - m_DataPoints(m, k));
        workValue = std::cos(2.0 * M_PI * scalarProduct);
      }
      else
      {
        workValue = 1.0;
        for (unsigned int k = 0;k < m_DomainDimension;++k)
        {
          if (std::abs(kVector[k]) < m_Epsilon)
            continue;
          workValue *= std::cos(2.0 * M_PI * kVector[k] * (m_DataPoints(l, k) - m_DataPoints(m, k)));
        }
      }

      workValue *= m_InternalLMatrix(m_PointLabels[l] - 1, m_PointLabels[m] - 1);

      m_DataLMatrix(l, m) += workValue;
      if (l != m)
        m_DataLMatrix(m, l) += workValue;
    }
  }
}

double BaseLogLikelihood::Evaluate(const arma::mat& x)
{
  this->SetModelParameters(x);

  if (m_FirstAmplitude < m_Epsilon || m_FirstAmplitude > 1.0)
    return(DBL_MAX);

  m_LogSpectrum = 0.0;
  m_DataLMatrix.fill(0.0);
  m_TraceValue = 0.0;
  m_ContinueLoop = true;
  double weightValue = 1.0;
  unsigned int counter = 0;

  if (m_SplitSummation)
  {
    for (unsigned int i = 0;i <= m_DomainDimension;++i)
    {
      weightValue = (double)std::pow(2, i);

      for (unsigned int j = 0;j < m_OptimizedIntegerGrid[i].size();++j)
      {
        ++counter;

        this->IncrementSummation(m_OptimizedIntegerGrid[i][j].first, m_OptimizedIntegerGrid[i][j].second, weightValue);

        if (!m_ContinueLoop)
          break;
      }

      if (!m_ContinueLoop)
        break;
    }
  }
  else
  {
    for (unsigned int i = 0;i < m_IntegerGrid.size();++i)
    {
      ++counter;

      this->IncrementSummation(m_IntegerGrid[i].first, m_IntegerGrid[i].second, weightValue);

      if (!m_ContinueLoop)
        break;
    }
  }

  // if (m_UseVerbose)
  Rcpp::Rcout << "* Number of k-vector used:      " << counter << std::endl;

  for (unsigned int i = 0;i < m_NumberOfPoints;++i)
  {
    for (unsigned int j = i;j < m_NumberOfPoints;++j)
    {
      if (std::abs(m_DataLMatrix(i, j)) < m_Epsilon)
      {
        m_DataLMatrix(i, j) = 0.0;
        if (i != j)
          m_DataLMatrix(j, i) = 0.0;
      }
    }
  }

  double sign;
  arma::log_det(m_LogDeterminant, sign, m_DataLMatrix);

  // if (m_UseVerbose)
  Rcpp::Rcout << "* Log-spectrum: " << m_LogSpectrum << std::endl;
  Rcpp::Rcout << "* Log-determinant: " << m_LogDeterminant << std::endl;

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

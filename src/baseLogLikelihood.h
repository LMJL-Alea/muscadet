#ifndef BASELOGLIKELIHOOD_H
#define BASELOGLIKELIHOOD_H

#include <RcppArmadillo.h>

class BaseLogLikelihood
{
public:
  BaseLogLikelihood()
  {
    m_NumberOfPoints = 1;
    m_NumberOfMarks = 1;
    m_DomainDimension = 1;
    m_MaximalNumberOfKVectors = 1;
    m_NumberOfKVectors = 1;
    m_TruncationIndex = 512;
    m_ActualTruncationIndex = 1;
    m_NumberOfThreads = 1;
    m_NumberOfParameters = 1;

    m_RelativeTolerance = 0.99;
    m_LogSpectrum = 0.0;
    m_LogDeterminant = 0.0;
    m_DomainVolume = 1.0;

    m_UseVerbose = false;
    m_UseFixedMarginalParameters = false;

    m_DataPoints.reset();
    m_WorkingEigenValues.reset();
    m_DeltaDiagonal.reset();
    m_LMatrixSum.reset();
    m_PointLabels.reset();
    m_DataLMatrix.reset();
    m_InternalLMatrix.reset();
    m_WorkingEigenVectors.reset();
    m_ListOfInternalLMatrices.reset();
  }

  ~BaseLogLikelihood() {}

  void SetInputData(
      const arma::mat &points,
      const arma::vec &lb,
      const arma::vec &ub,
      const arma::uvec &labels,
      const Rcpp::DataFrame &ndGrid,
      const unsigned int N,
      const Rcpp::Nullable<arma::vec> &marginal_parameters = R_NilValue
  );

  // Setter/getter for number of threads
  void SetNumberOfThreads(const unsigned int val) {m_NumberOfThreads = val;}
  unsigned int GetNumberOfThreads() const {return m_NumberOfThreads;}

  // Setter/getter for number of parameters
  void SetNumberOfParameters(const unsigned int val) {m_NumberOfParameters = val;}
  unsigned int GetNumberOfParameters() const {return m_NumberOfParameters;}

  // Setter/getter for alpha1
  void SetFirstAlpha(const double val);
  double GetFirstAlpha() const {return m_FirstAlpha;}

  // Setter/getter for alpha2
  void SetSecondAlpha(const double val);
  double GetSecondAlpha() const {return m_SecondAlpha;}

  // Setter/getter for alpha12
  void SetCrossParameters(const double alpha12, const double tau);
  double GetCrossAlpha() const {return m_CrossAlpha;}
  double GetCorrelation() const {return m_Correlation;}

  void SetUseVerbose(const bool &val) {m_UseVerbose = val;}
  bool GetUseVerbose() const {return m_UseVerbose;}

  // Return the objective function f(x) for the given x.
  double GetValue(const arma::vec& x);

  double GetFirstIntensity() const {return m_FirstIntensity;}
  double GetSecondIntensity() const {return m_SecondIntensity;}
  unsigned int GetDomainDimension() const {return m_DomainDimension;}

  double GetSquaredCrossAmplitudeUpperBound() const;
  virtual double GetCrossAlphaLowerBound() const = 0;
  virtual double RetrieveIntensityFromParameters(
      const double amplitude,
      const double alpha,
      const unsigned int dimension) = 0;
  virtual double RetrieveAlphaFromParameters(
      const double amplitude,
      const double intensity,
      const unsigned int dimension
  ) = 0;
  virtual double RetrieveAmplitudeFromParameters(
      const double intensity,
      const double alpha,
      const unsigned int dimension) = 0;

protected:

  //! Generic functions to be implemented in each child class
  virtual double GetK11Value(const double squaredNorm) = 0;
  virtual double GetK12Value(const double squaredNorm) = 0;
  virtual double GetK22Value(const double squaredNorm) = 0;
  double GetFirstAmplitude() {return m_FirstAmplitude;}
  double GetSecondAmplitude() {return m_SecondAmplitude;}
  double GetCrossAmplitude() {return m_CrossAmplitude;}

private:
  void SetModelParameters(const arma::vec &params);
  bool CheckModelParameters();
  void ComputeLogSpectrum();
  void ComputeLogDeterminant();

  unsigned int m_NumberOfPoints;
  unsigned int m_NumberOfMarks;
  unsigned int m_DomainDimension;
  unsigned int m_MaximalNumberOfKVectors;
  unsigned int m_NumberOfKVectors;
  unsigned int m_TruncationIndex;
  unsigned int m_ActualTruncationIndex;
  unsigned int m_NumberOfThreads;
  unsigned int m_NumberOfParameters;

  double m_RelativeTolerance;
  double m_LogSpectrum;
  double m_LogDeterminant;
  double m_DomainVolume;
  double m_FirstAlpha;
  double m_SecondAlpha;
  double m_CrossAlpha;
  double m_FirstIntensity;
  double m_SecondIntensity;
  double m_Correlation;
  double m_FirstAmplitude;
  double m_CrossAmplitude;
  double m_SecondAmplitude;

  bool m_UseVerbose;
  bool m_UseFixedMarginalParameters;

  arma::vec m_WorkingEigenValues;
  arma::vec m_DeltaDiagonal;

  arma::uvec m_PointLabels;

  arma::mat m_DataPoints;
  arma::mat m_DataLMatrix;
  arma::mat m_InternalLMatrix;
  arma::mat m_WorkingEigenVectors;
  arma::mat m_LMatrixSum;

  arma::cube m_ListOfInternalLMatrices;

  Rcpp::IntegerMatrix m_KGrid;
  Rcpp::NumericVector m_KSquaredNorms;
  Rcpp::IntegerVector m_KWeights;

  static const double m_Epsilon;
};

#endif /* BASELOGLIKELIHOOD_H */

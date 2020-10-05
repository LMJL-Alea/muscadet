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
    m_ActualTruncationIndex = 1;
    m_TruncationIndex = 512;

    m_TraceValue = 0.0;
    m_RelativeTolerance = 0.99;
    m_LogSpectrum = 0.0;
    m_LogDeterminant = 0.0;
    m_DomainVolume = 1.0;

    m_UseVerbose = false;

    m_ParameterLowerBounds.reset();
    m_ParameterUpperBounds.reset();
    m_WorkingEigenValues.reset();
    m_GradientIntegral.reset();
    m_GradientLogDeterminant.reset();
    m_ConstraintVector.reset();
    m_DeltaDiagonal.reset();
    m_LMatrixSum.reset();
    m_PointLabels.reset();
    m_DataLMatrix.reset();
    m_InternalLMatrix.reset();
    m_WorkingEigenVectors.reset();
    m_DataPoints.reset();
    m_ListOfInternalLMatrices.reset();
    m_CosineMatrix.reset();
    m_CosineValues.reset();
  }

  ~BaseLogLikelihood() {}

  void SetInputData(
      const arma::mat &points,
      const arma::vec &lb,
      const arma::vec &ub,
      const arma::uvec &labels,
      const Rcpp::DataFrame &ndGrid,
      const unsigned int N
  );

  // to be used in optimizer class
  unsigned int GetNumberOfParameters();
  arma::mat GetInitialPoint();
  arma::vec GetParameterLowerBounds() {return m_ParameterLowerBounds;}
  arma::vec GetParameterUpperBounds() {return m_ParameterUpperBounds;}

  // Setter/getter for alpha1
  void SetFirstAlpha(const double val);
  double GetFirstAlpha() {return m_FirstAlpha;}

  // Setter/getter for alpha2
  void SetSecondAlpha(const double val);
  double GetSecondAlpha() {return m_SecondAlpha;}

  // Setter/getter for alpha12
  void SetCrossAlpha(const double val) {m_CrossAlpha = val;}
  double GetCrossAlpha() {return m_CrossAlpha;}

  // Setter/getter for tau
  void SetCorrelation(const double val) {m_Correlation = val;}
  double GetCorrelation() {return m_Correlation;}

  void SetUseVerbose(const bool &val) {m_UseVerbose = val;}

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

  // Return the objective function f(x) for the given x.
  double Evaluate(const arma::mat& x);

  // Compute the gradient of f(x) for the given x and store the result in g.
  void Gradient(const arma::mat& x, arma::mat& g);

  // OPTIONAL: this may be implemented in addition to---or instead
  // of---Evaluate() and Gradient().  If this is the only function implemented,
  // implementations of Evaluate() and Gradient() will be automatically
  // generated using template metaprogramming.  Often, implementing
  // EvaluateWithGradient() can result in more efficient optimizations.
  //
  // Given parameters x and a matrix g, return the value of f(x) and store
  // f'(x) in the provided matrix g.  g should have the same size (rows,
  // columns) as x.
  double EvaluateWithGradient(const arma::mat& x, arma::mat& g);

  // Get the number of constraints on the objective function.
  size_t NumConstraints() {return 5;}

  // Evaluate constraint i at the parameters x.  If the constraint is
  // unsatisfied, DBL_MAX should be returned.  If the constraint is satisfied,
  // any real value can be returned.  The optimizer will add this value to its
  // overall objective that it is trying to minimize.  (So, a hard constraint
  // can just return 0 if it's satisfied.)
  double EvaluateConstraint(const size_t i, const arma::mat& x);

  // Evaluate the gradient of constraint i at the parameters x, storing the
  // result in the given matrix g.  If this is a hard constraint you can set
  // the gradient to 0.  If the constraint is not satisfied, it could be
  // helpful to set the gradient in such a way that the gradient points in the
  // direction where the constraint would be satisfied.
  void GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g);

protected:
  double GetCrossAmplitudeNormalizationFactor();
  //! Generic functions to be implemented in each child class
  virtual double GetCrossAlphaLowerBound() = 0;
  virtual double GetK11Value(const double squaredNorm) = 0;
  virtual double GetK12Value(const double squaredNorm) = 0;
  virtual double GetK22Value(const double squaredNorm) = 0;
  double GetFirstAmplitude() {return m_FirstAmplitude;}
  double GetSecondAmplitude() {return m_SecondAmplitude;}
  double GetCrossAmplitude() {return m_CrossAmplitude;}
  unsigned int GetDomainDimension() {return m_DomainDimension;}

private:
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters();
  void ComputeAll();
  void ComputeLogSpectrum();
  void ComputeLogDeterminant();

  unsigned int m_NumberOfPoints;
  unsigned int m_NumberOfMarks;
  unsigned int m_DomainDimension;
  unsigned int m_MaximalNumberOfKVectors;
  unsigned int m_NumberOfKVectors;
  unsigned int m_ActualTruncationIndex;

  int m_TruncationIndex;

  double m_TraceValue;
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
  double m_CrossBeta;
  double m_FirstAmplitude;
  double m_CrossAmplitude;
  double m_SecondAmplitude;
  double m_NormalizedCrossAmplitude;
  double m_InverseCrossAlpha;

  bool m_UseVerbose;

  arma::vec m_ParameterLowerBounds;
  arma::vec m_ParameterUpperBounds;
  arma::vec m_WorkingEigenValues;
  arma::vec m_GradientIntegral;
  arma::vec m_GradientLogDeterminant;
  arma::vec m_ConstraintVector;
  arma::vec m_DeltaDiagonal;

  arma::uvec m_PointLabels;

  arma::mat m_DataLMatrix;
  arma::mat m_InternalLMatrix;
  arma::mat m_WorkingEigenVectors;
  arma::mat m_DataPoints;
  arma::mat m_CosineValues;
  arma::mat m_LMatrixSum;
  arma::mat m_CosineMatrix;

  arma::cube m_ListOfInternalLMatrices;

  Rcpp::IntegerMatrix m_KGrid;
  Rcpp::NumericVector m_KSquaredNorms;
  Rcpp::IntegerVector m_KWeights;

  static const double m_Epsilon;
};

#endif /* BASELOGLIKELIHOOD_H */

#ifndef BASELOGLIKELIHOOD_H
#define BASELOGLIKELIHOOD_H

#include <RcppArmadillo.h>

class BaseLogLikelihood
{
public:
  using KVectorType = std::vector<double>;
  using KVectorPairType = std::pair<double,KVectorType>;

  BaseLogLikelihood()
  {
    m_TruncationIndex = 50;
    m_RelativeTolerance = 0.99;
    m_SplitSummation = true;
  }

  ~BaseLogLikelihood() {}

  void SetInputData(
      const arma::mat &points,
      const arma::vec &lb,
      const arma::vec &ub,
      const arma::uvec &labels
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
  void SetCrossAlpha(const double val);
  double GetCrossAlpha() {return m_CrossAlpha;}

  // Setter/getter for tau
  void SetCorrelation(const double val);
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

  void SetTruncationIndex(const int index) {m_TruncationIndex = index;}

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
  void GenerateCombinations(
      const unsigned int N,
      const unsigned int K,
      std::vector<std::vector<unsigned int> > &resVector
  );
  void SetIntegerGrid();
  void IncrementSummation(
      const double kSquaredNorm,
      const KVectorType &kVector,
      const double weight
  );

  unsigned int m_NumberOfPoints;

  arma::vec m_ParameterLowerBounds, m_ParameterUpperBounds;

  arma::mat m_DataLMatrix;
  arma::mat m_InternalLMatrix;
  double m_TraceValue;
  bool m_ContinueLoop;
  arma::vec m_WorkingEigenValues;
  arma::mat m_WorkingEigenVectors;
  double m_RelativeTolerance;
  bool m_SplitSummation;

  std::vector<KVectorPairType> m_IntegerGrid;
  std::vector<std::vector<KVectorPairType> > m_OptimizedIntegerGrid;

  double m_LogSpectrum, m_LogDeterminant;
  arma::vec m_GradientIntegral, m_GradientLogDeterminant;
  arma::uvec m_PointLabels;
  arma::vec m_ConstraintVector;
  double m_DomainVolume;
  arma::vec m_DeltaDiagonal;
  int m_TruncationIndex;
  arma::mat m_DataPoints;
  unsigned int m_NumberOfMarks;

  //! Generic variables used by all models and needed in each child class
  unsigned int m_DomainDimension;

  // Original parameters
  double m_FirstAlpha, m_SecondAlpha, m_CrossAlpha;
  double m_FirstIntensity, m_SecondIntensity, m_Correlation;

  // Variables for optimized parameters
  double m_CrossBeta;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
  double m_NormalizedCrossAmplitude;
  double m_InverseCrossAlpha;

  bool m_UseVerbose;

  static const double m_Epsilon;
};

#endif /* BASELOGLIKELIHOOD_H */

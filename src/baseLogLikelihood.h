#include "integrandFunctions.h"
#include <RcppEnsmallen.h>

class BaseLogLikelihood
{
public:
  typedef std::vector  <std::vector <int> > NeighborhoodType;
  typedef BaseIntegrand::KFunctionType KFunctionType;

  BaseLogLikelihood()
  {
    m_FirstAlpha = NA_REAL;
    m_SecondAlpha = NA_REAL;
    m_CrossAlpha = NA_REAL;
    m_FirstIntensity = NA_REAL;
    m_SecondIntensity = NA_REAL;
    m_CrossIntensity = NA_REAL;
    m_FirstAmplitude = NA_REAL;
    m_SecondAmplitude = NA_REAL;
    m_CrossAmplitude = NA_REAL;
    m_EstimateFirstAmplitude = true;
    m_EstimateSecondAmplitude = true;
    m_EstimateCrossAmplitude = true;
    m_EstimateFirstAlpha = true;
    m_EstimateSecondAlpha = true;
    m_EstimateCrossAlpha = true;

    m_DomainDimension = 1;
    m_DomainVolume = 1.0;
    m_UsePeriodicDomain = true;
    m_Modified = true;
    m_Integral = 0.0;
    m_LogDeterminant = 0.0;
  }

  ~BaseLogLikelihood() {}

  void SetInputs(
      const arma::mat &points,
      const arma::uvec &labels,
      const arma::vec &lb,
      const arma::vec &ub
  );
  void SetUsePeriodicDomain(const bool x) {m_UsePeriodicDomain = x;}
  arma::mat GetInitialPoint();

  void SetFirstAlpha(const double x);
  void SetSecondAlpha(const double x);
  void SetCrossAlpha(const double x);
  void SetFirstAmplitude(const double x);
  void SetSecondAmplitude(const double x);
  void SetCrossAmplitude(const double x);

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
  //! Generic functions to be implemented in each child class
  virtual double EvaluateLFunction(
      const double sqDist,
      const double amplitude,
      const double amplitude12,
      const double alpha,
      const double alpha12,
      const double l12Value,
      const unsigned int dimension) = 0;
  virtual double EvaluateL12Function(
      const double sqDist,
      const double amplitude1,
      const double amplitude2,
      const double amplitude12,
      const double alpha12,
      const unsigned int dimension) = 0;
  virtual double RetrieveIntensityFromParameters(
      const double amplitude,
      const double alpha,
      const unsigned int dimension) = 0;
  virtual double RetrieveAlphaFromParameters(
      const double amplitude,
      const double intensity,
      const unsigned int dimension
  ) = 0;
  virtual bool EvaluateAlphaConstraint(
      const double firstAlpha,
      const double secondAlpha,
      const double crossAlpha) = 0;
  virtual KFunctionType GetKFunction() = 0;
  double GetBesselJRatio(const double sqDist, const double alpha, const unsigned int dimension);

private:
  //! Helper functions for periodizing the domain
  unsigned int GetNumberOfParameters();
  void SetNeighborhood(const unsigned int n);
  std::vector<arma::rowvec> GetTrialVectors(const arma::rowvec &x, const arma::vec &lb, const arma::vec &ub);
  void SetModelParameters(const arma::mat &params);
  bool CheckModelParameters();
  double GetIntegral();
  double GetLogDeterminant();

  //! Generic variables used by all models but not needed in child classes
  double m_Integral, m_LogDeterminant;
  arma::vec m_GradientIntegral, m_GradientLogDeterminant;
  NeighborhoodType m_Neighborhood;
  bool m_UsePeriodicDomain;
  unsigned int m_SampleSize;
  arma::mat m_DistanceMatrix;
  arma::uvec m_PointLabels;
  arma::vec m_ConstraintVector;
  bool m_Modified;
  double m_DomainVolume;

  //! Generic variables used by all models and needed in each child class
  unsigned int m_DomainDimension;
  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
  bool m_EstimateFirstAmplitude, m_EstimateCrossAmplitude, m_EstimateSecondAmplitude;
  bool m_EstimateFirstAlpha, m_EstimateCrossAlpha, m_EstimateSecondAlpha;

  static const double m_Epsilon;
};

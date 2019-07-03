#pragma once

#include <RcppEnsmallen.h>

class BaseLogLikelihood
{
public:
  typedef std::vector  <std::vector <int> > NeighborhoodType;

  BaseLogLikelihood()
  {
    m_DomainDimension = 1;
    m_DomainVolume = 1.0;
    m_UsePeriodicDomain = true;
    m_Modified = true;
    m_FirstAlpha = 0.0;
    m_CrossAlpha = 0.0;
    m_SecondAlpha = 0.0;
    m_FirstIntensity = 0.0;
    m_CrossIntensity = 0.0;
    m_SecondIntensity = 0.0;
    m_FirstAmplitude = 0.0;
    m_CrossAmplitude = 0.0;
    m_SecondAmplitude = 0.0;
    m_Integral = 0.0;
    m_LogDeterminant = 0.0;
    m_GradientIntegral.set_size(4);
    m_GradientLogDeterminant.set_size(4);
  }

  ~BaseLogLikelihood() {}

  void SetInputs(
      const arma::mat &points,
      const arma::vec &labels,
      const arma::vec &lb,
      const arma::vec &ub,
      const double rho1,
      const double rho2
  );
  void SetUsePeriodicDomain(const bool x) {m_UsePeriodicDomain = x;}

  // Return the objective function f(x) for the given x.
  double Evaluate(const arma::mat& x);

  // Compute the gradient of f(x) for the given x and store the result in g.
  void Gradient(const arma::mat& x, arma::mat& g);

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

  arma::mat GetDistanceMatrix() {return m_DistanceMatrix;}

protected:
  void SetModelParameters(const arma::mat &params);
  virtual bool CheckModelParameters() = 0;
  virtual double GetIntegral() = 0;
  virtual double GetLogDeterminant() = 0;

  double m_FirstAlpha, m_CrossAlpha, m_SecondAlpha;
  double m_FirstIntensity, m_CrossIntensity, m_SecondIntensity;
  double m_FirstAmplitude, m_CrossAmplitude, m_SecondAmplitude;
  arma::vec m_GradientIntegral, m_GradientLogDeterminant;
  unsigned int m_DomainDimension;
  unsigned int m_SampleSize;
  arma::mat m_DistanceMatrix;
  arma::vec m_PointLabels;
  arma::vec m_ConstraintVector;

private:
  void SetNeighborhood(const unsigned int n);
  std::vector<arma::rowvec> GetTrialVectors(const arma::rowvec &x, const arma::vec &lb, const arma::vec &ub);

  double m_DomainVolume;
  bool m_Modified;
  double m_Integral, m_LogDeterminant;
  arma::vec m_DomainLowerBounds, m_DomainUpperBounds;
  NeighborhoodType m_Neighborhood;
  bool m_UsePeriodicDomain;
};

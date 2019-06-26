#include <RcppEnsmallen.h>

class BaseLogLikelihood
{
public:
  BaseLogLikelihood()
  {
    m_DataVolume = 1.0;
    m_Modified = true;
    m_Alpha1 = 0.0;
    m_Alpha12 = 0.0;
    m_Alpha2 = 0.0;
    m_Covariance = 0.0;
    m_Integral = 0.0;
    m_LogDeterminant = 0.0;
    m_GradientIntegral.set_size(4);
    m_GradientLogDeterminant.set_size(4);
  }

  ~BaseLogLikelihood() {}

  void SetInputs(const arma::mat &points, const double volume = 1.0);

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
  virtual double EvaluateConstraint(const size_t i, const arma::mat& x) = 0;

  // Evaluate the gradient of constraint i at the parameters x, storing the
  // result in the given matrix g.  If this is a hard constraint you can set
  // the gradient to 0.  If the constraint is not satisfied, it could be
  // helpful to set the gradient in such a way that the gradient points in the
  // direction where the constraint would be satisfied.
  void GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g)
  {
    g.set_size(x.n_cols, 1);
    g.fill(0.0);
  }

  arma::mat GetDistanceMatrix() {return m_DistanceMatrix;}

protected:
  void SetModelParameters(const arma::mat &params);
  virtual double GetIntegral() = 0;
  virtual double GetLogDeterminant() = 0;

  double m_Alpha1, m_Alpha12, m_Alpha2, m_Covariance;
  double m_Intensity1, m_Intensity2;
  arma::vec m_GradientIntegral, m_GradientLogDeterminant;
  unsigned int m_DataDimension;
  unsigned int m_SampleSize;
  arma::mat m_DistanceMatrix;
  arma::vec m_PointLabels;

private:
  double m_DataVolume;
  bool m_Modified;
  double m_Integral, m_LogDeterminant;
};

class GaussianLogLikelihood : public BaseLogLikelihood
{
public:
  double EvaluateConstraint(const size_t i, const arma::mat& x);

protected:
  double GetIntegral();
  double GetLogDeterminant();
};

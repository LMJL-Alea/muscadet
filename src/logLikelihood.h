#include <RcppEnsmallen.h>

class BaseLogLikelihood
{
public:
  BaseLogLikelihood()  {}
  ~BaseLogLikelihood() {}

  void SetInputs(const arma::mat &points, const double volume);
  double GetIntensity1() {return m_Intensity1;}
  double GetIntensity2() {return m_Intensity2;}
  arma::mat GetDistanceMatrix() {return m_DistanceMatrix;}
  arma::vec GetPointLabels() {return m_PointLabels;}
  double GetDataVolume() {return m_DataVolume;}
  unsigned int GetDataDimension() {return m_DataDimension;}
  unsigned int GetSampleSize() {return m_SampleSize;}

  virtual double GetNormalizationFactor(const arma::mat &params) = 0;
  virtual double GetLogDeterminant(const arma::mat &params, arma::vec &grad) = 0;

  // Return the objective function f(x) for the given x.
  double Evaluate(const arma::mat& x);

  // Compute the gradient of f(x) for the given x and store the result in g.
  void Gradient(const arma::mat& x, arma::mat& g);

  // Get the number of constraints on the objective function.
  unsigned int NumConstraints() {return 5;}

  // Evaluate constraint i at the parameters x.  If the constraint is
  // unsatisfied, DBL_MAX should be returned.  If the constraint is satisfied,
  // any real value can be returned.  The optimizer will add this value to its
  // overall objective that it is trying to minimize.  (So, a hard constraint
  // can just return 0 if it's satisfied.)
  virtual double EvaluateConstraint(const unsigned int i, const arma::mat& x) = 0;

  // Evaluate the gradient of constraint i at the parameters x, storing the
  // result in the given matrix g.  If this is a hard constraint you can set
  // the gradient to 0.  If the constraint is not satisfied, it could be
  // helpful to set the gradient in such a way that the gradient points in the
  // direction where the constraint would be satisfied.
  void GradientConstraint(const unsigned i, const arma::mat& x, arma::mat& g) {g.fill(0.0);}

protected:
  double m_Intensity1, m_Intensity2;
  unsigned int m_DataDimension;
  unsigned int m_SampleSize;
  arma::mat m_DistanceMatrix;
  arma::vec m_PointLabels;

private:
  double m_DataVolume;
};

class GaussianLogLikelihood : public BaseLogLikelihood
{
public:
  double GetNormalizationFactor(const arma::mat &params);
  double GetLogDeterminant(const arma::mat &params, arma::vec &grad);

  double EvaluateConstraint(const unsigned int i, const arma::mat& x);
};

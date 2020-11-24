#include <RcppArmadillo.h>
#include "baseLogLikelihood.h"
#include "baseOptimizerClass.h"

class DeterminantalPointProcess
{
public:
  DeterminantalPointProcess() {}
  ~DeterminantalPointProcess() {}
  void SetLikelihoodModel(const std::string &val);
  void SetOptimizer(const std::string &val);
  Rcpp::List Fit(
      const arma::mat &points,
      const arma::vec &lower_bound,
      const arma::vec &upper_bound,
      const Rcpp::DataFrame &nd_grid,
      const Rcpp::Nullable<arma::uvec> &marks = R_NilValue,
      const unsigned int num_threads = 1,
      const unsigned int N = 50,
      const bool use_verbose = false
  ) const;

private:
  std::shared_ptr<BaseLogLikelihood> m_LikelihoodPointer;
  std::shared_ptr<BaseOptimizerFunction> m_OptimizerPointer;
};

RCPP_MODULE(DPP) {
  using namespace Rcpp;
  class_<DeterminantalPointProcess>("DeterminantalPointProcess")
    .default_constructor("Default constructor")
    .method("SetLikelihoodModel", &DeterminantalPointProcess::SetLikelihoodModel)
    .method("SetOptimizer",       &DeterminantalPointProcess::SetOptimizer)
    .method("Fit",                &DeterminantalPointProcess::Fit)
  ;
}

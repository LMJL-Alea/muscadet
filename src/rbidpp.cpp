#include "rbidpp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

arma::mat22 GetGaussianKernel(const double rSq,
                              const double rho1,
                              const double rho2,
                              const double alpha1Sq,
                              const double alpha2Sq,
                              const double alpha12Sq,
                              const double tau)
{
  arma::mat22 kernelMatrix;
  double piVal = M_PI;
  double rho12 = tau * std::sqrt(rho1* rho2);
  kernelMatrix(0, 0) = rho1 * alpha1Sq * piVal * std::exp(-piVal * piVal * alpha1Sq * rSq);
  kernelMatrix(0, 1) = rho12 * alpha12Sq * piVal * std::exp(-piVal * piVal * alpha12Sq * rSq);
  kernelMatrix(1, 0) = kernelMatrix(0, 1);
  kernelMatrix(1, 1) = rho2 * alpha2Sq * piVal * std::exp(-piVal * piVal * alpha2Sq * rSq);
  return kernelMatrix;
}

Rcpp::List rbidpp_impl(const int N,
                       const double L,
                       const double rho1,
                       const double rho2,
                       const double alpha1,
                       const double alpha2,
                       const double alpha12,
                       const double tau,
                       const unsigned int nbThreads)
{
  std::vector<arma::imat> threadedKIndices(nbThreads);
  std::vector<arma::mat> threadedEigenVectors(nbThreads);
  for (unsigned int i = 0;i < nbThreads;++i)
  {
    threadedKIndices[i] = arma::imat(0, 2);
    threadedEigenVectors[i] = arma::mat(2, 0);
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nbThreads)
#endif

  for (int i = -N;i <= N;++i)
  {
    arma::vec2 workEigenValues;
    arma::mat22 workEigenVectors;
    arma::mat22 workKernelMatrix;
    Rcpp::NumericVector workBernoulli(1);
    arma::irowvec2 workIndices;
    unsigned int threadId = omp_get_thread_num();
    // printf("Hello world from omp thread %d\n", threadId);

    workIndices(0) = i;

    for (int j = -N;j <= N;++j)
    {
      workIndices(1) = j;

      workKernelMatrix = GetGaussianKernel(
        (i * i + j * j) / (L * L),
        rho1, rho2,
        alpha1 * alpha1,
        alpha2 * alpha2,
        alpha12 * alpha12,
        tau
      );

      arma::eig_sym(workEigenValues, workEigenVectors, workKernelMatrix);

      for (unsigned int k = 0;k < 2;++k)
      {
        workBernoulli = Rcpp::rbinom(1, 1, workEigenValues(k));

        if (workBernoulli(0) == 1)
        {
          threadedKIndices[threadId].insert_rows(threadedKIndices[threadId].n_rows, workIndices);
          threadedEigenVectors[threadId].insert_cols(threadedEigenVectors[threadId].n_cols, workEigenVectors.col(k));
        }
      }
    }
  }

  // Rcpp::Rcout<< "Done with parallel loop" <<std::endl;

  arma::imat kIndices(0, 2);
  arma::mat eigenVectors(2, 0);

  for (unsigned int i = 0;i < nbThreads;++i)
  {
    kIndices.insert_rows(kIndices.n_rows, threadedKIndices[i]);
    eigenVectors.insert_cols(eigenVectors.n_cols, threadedEigenVectors[i]);
  }

  // Rcpp::Rcout<< "Done aggregating threaded data" <<std::endl;

  return Rcpp::List::create(
    Rcpp::Named("kkindex") = kIndices,
    Rcpp::Named("V") = eigenVectors
  );
}

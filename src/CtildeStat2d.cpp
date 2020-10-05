#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix CtildeStat2d_cpp(NumericVector x,
                               NumericVector y,
                               double lambdao,
                               NumericVector lambdaa,
                               NumericVector lambda,
                               IntegerVector k1a,
                               IntegerVector k2a )
{
  unsigned int nC = x.size();
  unsigned int na = k1a.size();
  unsigned int nn = na / 2;
  NumericVector cos1(nn + 1);
  NumericVector cos2(nn + 1);
  NumericMatrix C(nC, nC);

  double diagsum = lambdao + 2 * sum(lambdaa) + 4 * sum(lambda);

  for (unsigned int ii = 0;ii < nC - 1;++ii)
  {
    // Assign the diagonal
    C(ii, ii) = diagsum;

    // Assign upper triangle (and copy it to lower)
    for (unsigned int jj = ii + 1;jj < nC;++jj)
    {
      double xdiff = 2.0 * M_PI * (x[ii] - x[jj]);
      double ydiff = 2.0 * M_PI * (y[ii] - y[jj]);

      double sum1 = 0.0;
      for (unsigned int kk = 0;kk < na;++kk)
      	sum1 += lambdaa[kk] * std::cos(k1a[kk] * xdiff + k2a[kk] * ydiff);

      cos1[0] = 1.0;
      cos2[0] = 1.0;
      cos1[1] = std::cos(xdiff);
      cos2[1] = std::cos(ydiff);
      for (unsigned int kk = 2;kk <= nn;++kk)
      {
        cos1[kk] = 2.0 * cos1[1] * cos1[kk - 1] - cos1[kk - 2];
        cos2[kk] = 2.0 * cos2[1] * cos2[kk - 1] - cos2[kk - 2];
      }

      double sum2 = 0.0;
      for (unsigned int kk = 0;kk < nn;++kk)
        for (unsigned int ll = 0;ll < nn;++ll)
          sum2 += lambda[ll + nn * kk] * cos1[ll + 1] * cos2[kk + 1];

      C(ii, jj) = lambdao + 2.0 * sum1 + 4.0 * sum2;
      C(jj, ii) = lambdao + 2.0 * sum1 + 4.0 * sum2;
    }
  }

  // Assign final diagonal element
  C(nC - 1, nC - 1) = diagsum;

  return C;
}

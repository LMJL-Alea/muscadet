% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcf.R
\name{fit_via_pcf}
\alias{fit_via_pcf}
\title{Estimation of Two-Mark Planar DPPs via the Pair Correlation Function}
\usage{
fit_via_pcf(
  X,
  model = "Gauss",
  rmin_alpha = 2,
  rmin_alpha12 = 2,
  rmin_tau = 2,
  q = 1,
  p = 2,
  divisor_marginal = "d",
  divisor_cross = "d",
  method = "profiling",
  B = 0L,
  full_bootstrap = TRUE,
  params = NULL
)
}
\arguments{
\item{X}{An object of class \link[spatstat.geom:ppp]{spatstat.geom::ppp} specifying a planar point
pattern.}

\item{model}{A string specifying the model to be fitted. Choices are
\code{"Gauss"} or \code{"Bessel"}. Defaults to \code{"Gauss"}.}

\item{rmin_alpha}{The lower bound on distances that should be taken into
account for estimating marginal alpha parameters (default: index 2).}

\item{rmin_alpha12}{The lower bound on distances that should be taken into
account for estimating the crossing alpha parameter (default: index 2).}

\item{rmin_tau}{The lower bound on distances that should be taken into
account for estimating the correlation (default: index 2).}

\item{q}{Power for pointwise evaluations of the PCF. Defaults to \code{0.5}.}

\item{p}{Power for the distance between empirical and moodel-based PCF
values. Defaults to \code{2}.}

\item{divisor_marginal}{Choice of divisor in the estimation formula. Choices
are \code{"r"} or \code{"d"}. See Section Empirical estimation of the pair
correlation function for more details. Defaults to \code{"d"}.}

\item{divisor_cross}{Choice of divisor in the estimation formula. Choices are
\code{"r"} or \code{"d"}. See Section Empirical estimation of the pair correlation
function for more details. Defaults to \code{"d"}.}

\item{method}{A character string specifying the estimation method between
\code{"profiling"} and \code{"direct"}. Defaults to \code{"profiling"}.}

\item{B}{An integer value specifying the number of samples to be generated in
the bootstrap procedure to approximate the distribution of the test
statistics when testing for absence of correlation between marks. Defaults
to \code{0L} which does not perform the test at all.}

\item{full_bootstrap}{A boolean specifying whether marginal parameters should
be re-estimated when computing the bootstrapped distribution of the tau
statistic. Defaults to \code{TRUE}.}

\item{params}{A length-4 numeric vector specifying values for the marginal
parameters if known. The order needs to be \code{rho1}, \code{rho2}, \code{alpha1} and
\code{alpha2}. Defaults to \code{NULL}, in which case, they are estimated.}
}
\value{
A list with the estimated model parameters in the following order:
\code{rho1}, \code{rho2}, \code{alpha1}, \code{alpha2}, \code{k12}, \code{alpha12} and \code{tau}. Additional
information pertaining to the test for absence of correlation between marks
are returned in the list as well.
}
\description{
Estimation of Two-Mark Planar DPPs via the Pair Correlation Function
}
\section{Empirical estimation of the pair correlation function}{

The empirical PCF is computed as a kernel estimate of the PCF in which the
contribution from an interpoint distance \eqn{d_{ij}} to the estimate of
\eqn{g(r)} is divided:
\itemize{
\item either by \eqn{r} using optional argument \code{divisor = "r"} in the functions
\code{\link[spatstat.explore:pcf]{spatstat.explore::pcf()}} and \code{\link[spatstat.explore:pcfcross]{spatstat.explore::pcfcross()}};
\item or by \eqn{d_{ij}} using optional argument \code{divisor = "d"} in the functions
\code{\link[spatstat.explore:pcf]{spatstat.explore::pcf()}} and \code{\link[spatstat.explore:pcfcross]{spatstat.explore::pcfcross()}}; it is intended
to improve the bias of the estimator when \eqn{r} is close to zero.
}
}

\examples{
fit_via_pcf(sim_gauss0[[1]])
}

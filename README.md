
<!-- README.md is generated from README.Rmd. Please edit that file -->

# muscadet

<!-- badges: start -->

[![check-standard](https://github.com/LMJL-Alea/muscadet/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/LMJL-Alea/muscadet/actions/workflows/check-standard.yaml)
[![pkgdown](https://github.com/LMJL-Alea/muscadet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/LMJL-Alea/muscadet/actions/workflows/pkgdown.yaml)
[![Codecov test
coverage](https://codecov.io/gh/LMJL-Alea/muscadet/branch/master/graph/badge.svg)](https://app.codecov.io/gh/LMJL-Alea/muscadet?branch=master)
[![test-coverage](https://github.com/LMJL-Alea/muscadet/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/LMJL-Alea/muscadet/actions/workflows/test-coverage.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/muscadet)](https://CRAN.R-project.org/package=muscadet)
<!-- badges: end -->

The goal of **muscadet** is to provide ways of simulating and estimating
planar point patterns according to various two-mark DPP models. It
currently provides support for Gaussian and Bessel DPPs. Estimation can
be performed either by minimizing the pair correlation function contrast
or by maximizing the likelihood.

## Installation

You can install the development version of **muscadet** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LMJL-Alea/muscadet")
```

## Example

The typical way of using **muscadet** for studying two-mark planar DPPs
is to use the
[`TwoMarkDPP`](https://lmjl-alea.github.io/muscadet/reference/TwoMarkDPP.html)
class:

``` r
library(muscadet)
mod <- TwoMarkDPP$new(model = "Gauss")
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: ?
#> • Intensity of the 2nd marginal DPP: ?
#> • Repulsion rate of the 1st marginal DPP: ?
#> • Repulsion rate of the 2nd marginal DPP: ?
#> • Repulsion rate between marks: ?
#> • Correlation between marks: ?
#> • Window size: ?
```

As one can see, the model is initially empty, in the sense that no
parameter is provided. You can start setting them manually one at a
time:

``` r
mod$first_intensity <- 100
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: 100
#> • Intensity of the 2nd marginal DPP: ?
#> • Repulsion rate of the 1st marginal DPP: (0, 0.056)
#> • Repulsion rate of the 2nd marginal DPP: ?
#> • Repulsion rate between marks: ?
#> • Correlation between marks: ?
#> • Window size: ?
```

You can notice how the model gets updated and gives you intervals of
feasible values for some parameters as soon as it is able to compute
this information given the set of already known parameters.

``` r
mod$second_intensity <- 100
mod$first_repulsion_rate <- 0.03
mod$second_repulsion_rate <- 0.03
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: 100
#> • Intensity of the 2nd marginal DPP: 100
#> • Repulsion rate of the 1st marginal DPP: 0.03
#> • Repulsion rate of the 2nd marginal DPP: 0.03
#> • Repulsion rate between marks: [0.03, Inf)
#> • Correlation between marks: ?
#> • Window size: ?
```

``` r
mod$cross_repulsion_rate <- 0.035
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: 100
#> • Intensity of the 2nd marginal DPP: 100
#> • Repulsion rate of the 1st marginal DPP: 0.03
#> • Repulsion rate of the 2nd marginal DPP: 0.03
#> • Repulsion rate between marks: 0.035
#> • Correlation between marks: [0, 0.735)
#> • Window size: ?
```

``` r
mod$between_mark_correlation <- 0.5
mod$window_size <- 1
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: 100
#> • Intensity of the 2nd marginal DPP: 100
#> • Repulsion rate of the 1st marginal DPP: 0.03
#> • Repulsion rate of the 2nd marginal DPP: 0.03
#> • Repulsion rate between marks: 0.035
#> • Correlation between marks: 0.5
#> • Window size: 1
```

Now the model is fully specified and you can simulate a point pattern
from it via the `$random()` method:

``` r
progressr::handlers("rstudio")
withr::with_seed(1234, {
  progressr::with_progress({
    df <- mod$random(n = 1)
  })
})
```

Note that the computations are parallelized over the sample size using
the [**futureverse**](https://future.futureverse.org) framework.

If you have an observed two-mark planar point pattern, you can also use
it to fit the model which automatically adjusts the parameter set and
bounds. This is achieved via the `$fit()` method:

``` r
mod$fit(df[[1]])
mod
#> 
#> ── Two-mark planar Gauss determinantal point process ───────────────────────────
#> • Intensity of the 1st marginal DPP: 105
#> • Intensity of the 2nd marginal DPP: 102
#> • Repulsion rate of the 1st marginal DPP: 0.034
#> • Repulsion rate of the 2nd marginal DPP: 0.037
#> • Repulsion rate between marks: 0.036
#> • Correlation between marks: 0.71
#> • Window size: 1
```

If you want finer control over the estimation methods, you can use the
functions
[`fit_via_pcf()`](https://lmjl-alea.github.io/muscadet/reference/fit_via_pcf.html)
and
[`fit_via_mle()`](https://lmjl-alea.github.io/muscadet/reference/fit_via_mle.html)
from outside the model class.

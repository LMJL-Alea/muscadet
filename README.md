
<!-- README.md is generated from README.Rmd. Please edit that file -->

# muscadet

<!-- badges: start -->

[![check-standard](https://github.com/LMJL-Alea/muscadet/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/LMJL-Alea/muscadet/actions/workflows/check-standard.yaml)
[![pkgdown](https://github.com/LMJL-Alea/muscadet/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/LMJL-Alea/muscadet/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

The goal of **muscadet** is to provide a set of functions for simulating
and estimating planar point patterns according to various two-mark DPP
models. It currently provides support for Gaussian and Bessel DPPs.

## Installation

You can install the development version of muscadet from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LMJL-Alea/muscadet")
```

## Example

You can simulate a planar point pattern following a two-mark DPP using
the `rbidpp()` function:

``` r
library(muscadet)
withr::with_seed(1234, {
  progressr::with_progress({
    df <- rbidpp(
      n = 1, 
      rho1 = 100, rho2 = 100, 
      alpha1 = 0.03, alpha2 = 0.03, 
      alpha12 = 0.035, tau = 0.5, 
      L = 1, model = "Gauss"
    )  
  })
})
```

You can estimate the parameters of a two-mark DPP from an observed
planar point pattern using the `estimate()` function:

``` r
estimate(df[[1]], model = "Gauss")
#> $rho1
#> [1] 105
#> 
#> $rho2
#> [1] 102
#> 
#> $alpha1
#> [1] 0.03351226
#> 
#> $alpha2
#> [1] 0.03741944
#> 
#> $k12
#> [1] 0.2910802
#> 
#> $tau
#> [1] 0.7096292
#> 
#> $alpha12
#> [1] 0.03551962
#> 
#> $fmin
#> [1] 0.001426345
#> 
#> $stat_np_obs
#> [1] 0.005231958
#> 
#> $stat_p_obs
#> [1] 0.7096292
#> 
#> $null_np_distr
#> NULL
#> 
#> $null_p_distr
#> NULL
#> 
#> $np_reject
#> NULL
#> 
#> $p_reject
#> NULL
```

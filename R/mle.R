#' Maximum Likelihood Estimator of Stationary Bivariate DPP
#'
#' @param X An n x d matrix storing n observed points in R^d.
#' @param labels An n-dimensional integer vector storing the labels of each
#'   observed point.
#' @param lb A d-dimensional numeric vector storing the lower bounds of the
#'   spatial domain (default: \code{rep(-0.5, ncol(X))}).
#' @param ub A d-dimensional numeric vector storing the upper bounds of the
#'   spatial domain (default: \code{rep( 0.5, ncol(X))}).
#' @param a1 Value of the first amplitude. If set to 0 (default), the parameter
#'   is estimated.
#' @param alpha1 Value of the first alpha If set to 0 (default), the parameter
#'   is estimated.
#' @param a2 Value of the second amplitude. If set to 0 (default), the parameter
#'   is estimated.
#' @param alpha2 Value of the second alpha If set to 0 (default), the parameter
#'   is estimated.
#' @param a12 Value of the cross amplitude. If set to 0 (default), the parameter
#'   is estimated.
#' @param alpha12 Value of the cross alpha If set to 0 (default), the parameter
#'   is estimated.
#'
#' @return A list as output from \code{\link[stats]{optim}}.
#' @name mle-dpp
#'
#' @examples
#' all <- bessel_tauSq0p5[[1]]
#' X <- cbind(all$x, all$y)
#' labels <- all$marks
#' rho1 <- rho2 <- 100
#' rho12 <- sqrt(0.5) * sqrt(rho1 * rho2)
#' alpha1 <- alpha2 <- 0.03
#' alpha12 <- 0.03
#' d <- 2
#' a1 <- rho1 * (2 * pi * alpha1 * alpha1 / d)^(d / 2) * gamma(d / 2 + 1)
#' a2 <- rho2 * (2 * pi * alpha2 * alpha2 / d)^(d / 2) * gamma(d / 2 + 1)
#' a12 <- rho12 * (2 * pi * alpha12 * alpha12 / d)^(d / 2) * gamma(d / 2 + 1)
#' mle_dpp_bessel(
#'   X = X,
#'   labels = labels,
#'   a1 = a1,
#'   a2 = a2,
#'   a12 = a12,
#'   alpha12 = alpha12
#' )
NULL

#' @rdname mle-dpp
#' @export
mle_dpp_gauss <- function(X, labels,
                          lb = rep(-0.5, ncol(X)),
                          ub = rep( 0.5, ncol(X)),
                          a1 = NA, alpha1 = NA,
                          a2 = NA, alpha2 = NA,
                          a12 = NA, alpha12 = NA) {
  x0 <- InitializeGauss(X, labels, lb, ub, a1, a2, a12, alpha1, alpha2, alpha12)

  optim(
    par = x0, fn = EvaluateGauss, method = "Nelder-Mead",
    control = list(warn.1d.NelderMead = FALSE),
    X = X, labels = labels, lb = lb, ub = ub,
    amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
    alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12
  )
}

#' @rdname mle-dpp
#' @export
mle_dpp_bessel <- function(X, labels,
                           lb = rep(-0.5, ncol(X)),
                           ub = rep( 0.5, ncol(X)),
                           a1 = NA, alpha1 = NA,
                           a2 = NA, alpha2 = NA,
                           a12 = NA, alpha12 = NA) {
  x0 <- InitializeBessel(X, labels, lb, ub, a1, a2, a12, alpha1, alpha2, alpha12)

  optim(
    par = x0, fn = EvaluateBessel, method = "Nelder-Mead",
    control = list(warn.1d.NelderMead = FALSE),
    X = X, labels = labels, lb = lb, ub = ub,
    amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
    alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12
  )
}

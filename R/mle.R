#' Maximum Likelihood Estimator of Stationary Bivariate DPPs
#'
#' @param X An n x d matrix storing n observed points in R^d.
#' @param labels An n-dimensional integer vector storing the labels of each
#'   observed point.
#' @param lb A d-dimensional numeric vector storing the lower bounds of the
#'   spatial domain (default: \code{rep(-0.5, ncol(X))}).
#' @param ub A d-dimensional numeric vector storing the upper bounds of the
#'   spatial domain (default: \code{rep( 0.5, ncol(X))}).
#' @param rho1 Value of the first intensity. If set, it requires toset also the
#'   \code{alpha1} parameter.
#' @param alpha1 Value of the first alpha If set to \code{NA} (default), the
#'   parameter is estimated.
#' @param rho2 Value of the second intensity. If set, it requires to set also
#'   the \code{alpha2} parameter.
#' @param alpha2 Value of the second alpha If set to \code{NA} (default), the
#'   parameter is estimated.
#' @param estimate_alpha A boolean specifying whether the marginal alpha's
#'   should be estimated (default: \code{TRUE}).
#'
#' @return A list as output from \code{\link[stats]{optim}}.
#' @name mle-dpp
#'
#' @examples
#' dpp <- bessel_tauSq0p5[[1]]
#' X <- cbind(dpp$x, dpp$y)
#' labels <- dpp$marks
#' rho1 <- rho2 <- 100
#' rho12 <- sqrt(0.5) * sqrt(rho1 * rho2)
#' alpha1 <- alpha2 <- 0.03
#' alpha12 <- 0.03
#' d <- 2
#' mle_dpp_bessel(
#'   X = X,
#'   labels = labels,
#'   rho1 = rho1,
#'   rho2 = rho2,
#'   alpha1 = alpha1,
#'   alpha2 = alpha2,
#'   estimate_alpha = FALSE
#' )
NULL

#' @rdname mle-dpp
#' @export
mle_dpp_gauss <- function(X, labels,
                          lb = rep(-0.5, ncol(X)),
                          ub = rep( 0.5, ncol(X)),
                          rho1 = NA, alpha1 = NA,
                          rho2 = NA, alpha2 = NA,
                          estimate_alpha = TRUE) {
  x0 <- InitializeGauss(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha)

  optim(
    par = x0, fn = EvaluateGauss, method = "Nelder-Mead",
    control = list(warn.1d.NelderMead = FALSE),
    X = X, labels = labels, lb = lb, ub = ub,
    rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2, estimate_alpha = estimate_alpha
  )
}

#' @rdname mle-dpp
#' @export
mle_dpp_bessel <- function(X, labels,
                           lb = rep(-0.5, ncol(X)),
                           ub = rep( 0.5, ncol(X)),
                           rho1 = NA, alpha1 = NA,
                           rho2 = NA, alpha2 = NA,
                           estimate_alpha = TRUE) {
  x0 <- InitializeBessel(X, labels, lb, ub, rho1, rho2, alpha1, alpha2, estimate_alpha)

  optim(
    par = x0, fn = EvaluateBessel, method = "Nelder-Mead",
    control = list(warn.1d.NelderMead = FALSE),
    X = X, labels = labels, lb = lb, ub = ub,
    rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2, estimate_alpha = estimate_alpha
  )
}

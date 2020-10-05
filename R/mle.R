#' Title
#'
#' @param theta
#' @param points
#' @param lower_bound
#' @param upper_bound
#' @param marks
#' @param N
#' @param use_verbose
#'
#' @return
#' @export
#'
#' @examples
my_mle <- function(theta, points, lower_bound, upper_bound, marks = NULL, N = 512, use_verbose = FALSE) {
  nd_grid <- generate_nd_grid(N, dim(points)[2])

  log_likelihood(
    theta = theta,
    points = points,
    lower_bound = lb,
    upper_bound = ub,
    nd_grid = nd_grid,
    marks = marks,
    N = N,
    use_verbose = use_verbose
  )
}

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
#' dpp <- sim[[1]]
#' rho1 <- rho2 <- 100
#' alpha1 <- alpha2 <- 0.03
#' alpha12 <- 0.05
#' tau <- 0.2
#' rho12 <- tau * sqrt(rho1 * rho2)
#' d <- 2
#' mle_dpp_bessel(
#'   X = X,
#'   lb = rep(0, ncol(X)),
#'   ub = rep(1, ncol(X))
#' )
#'
#' res <- lapply(sim, mle_dpp_bessel, lb = rep(0, ncol(X)), ub = rep(1, ncol(X)))
#' temp <- unlist(res)
#' par(mfrow = c(2, 2))
#' boxplot(temp[seq(2, 600, 6)], main = 'alpha1')
#' abline(h = 0.03)
#' boxplot(temp[seq(4, 600, 6)], main = 'alpha2')
#' abline(h = 0.03)
#' boxplot(temp[seq(5, 600, 6)], main = 'tau')
#' abline(h = 0.2)
#' boxplot(temp[seq(6, 600, 6)], main = 'alpha12')
#' abline(h = 0.05)
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

get_alpha <- function(k, rho, d) {
  (k / (rho * gamma(1 + d / 2)))^(1 / d) / sqrt(2 * pi / d)
}

get_k12 <- function(k12norm, k1, k2) {
  ub <- sqrt(max(min(k1 * k2, (1 - k1) * (1 - k2)), 0))
  k12norm * ub
}

get_alpha12bis <- function(beta12, alpha1, alpha2) {
  max(alpha1, alpha2) / beta12
}

get_tau <- function(k12, alpha12, rho1, rho2, d) {
  rho12 <- get_rho(k12, alpha12, d)
  rho12 /sqrt(rho1 * rho2)
}

#' @rdname mle-dpp
#' @export
mle_dpp_bessel <- function(X,
                           nlopt = "neldermead",
                           lb = rep(0, ncol(X)),
                           ub = rep(1, ncol(X)),
                           estimate_rho = TRUE,
                           init = NULL) {
  epsilon <- 1e-4

  labels <- X$marks
  if (is.null(init)) {
    rho2 <- spatstat::intensity(X)
    rho1 <- rho2[1]
    rho2 <- rho2[2]
  } else {
    rho1 <- init$rho1
    rho2 <- init$rho2
  }
  V <- spatstat::volume(X$window)
  X <- cbind(X$x, X$y)
  d <- ncol(X)

  if (!estimate_rho | is.null(init)) {
    # par is (k1, k2, k12norm, beta12)
    lbs <- c(
      epsilon,
      epsilon,
      0,
      epsilon
    )
    ubs <- c(
      1 - epsilon,
      1 - epsilon,
      1 - epsilon,
      1
    )

    # First, fit model with fixed rhos
    # Grab a good initial position
    if (is.null(init)) {
      fit <- nloptr::directL(
        fn = EvaluateBessel,
        lower = lbs,
        upper = ubs,
        X = X, labels = labels, lb = lb, ub = ub,
        rho1 = rho1, rho2 = rho2,
        original = TRUE
      )
      x0 <- fit$par
    } else {
      x0 <- c(
        init$alpha1,
        init$alpha2,
        get_k(init$tau * sqrt(rho1 * rho2), init$alpha12, d),
        init$tau
      )
    }

    if (nlopt == "bobyqa") {
      fit <- nloptr::bobyqa(
        x0 = x0,
        fn = EvaluateBessel,
        lower = lbs,
        upper = ubs,
        X = X, labels = labels, lb = lb, ub = ub,
        rho1 = rho1, rho2 = rho2
      )
    } else if (nlopt == "neldermead") {
      fit <- nloptr::neldermead(
        x0 = x0,
        fn = EvaluateBessel,
        lower = lbs,
        upper = ubs,
        X = X, labels = labels, lb = lb, ub = ub,
        rho1 = rho1, rho2 = rho2
      )
    } else if (nlopt == "sbplx") {
      fit <- nloptr::sbplx(
        x0 = x0,
        fn = EvaluateBessel,
        lower = lbs,
        upper = ubs,
        X = X, labels = labels, lb = lb, ub = ub,
        rho1 = rho1, rho2 = rho2
      )
    } else {
      # fit <- hydroPSO::hydroPSO(
      #   fn = "EvaluateBessel",
      #   lower = lbs,
      #   upper = ubs,
      #   X = X, labels = labels, lb = lb, ub = ub,
      #   rho1 = rho1, rho2 = rho2
      # )
      # fit <- GenSA::GenSA(
      #   lower = lbs,
      #   upper = ubs,
      #   fn = EvaluateBessel,
      #   X = X, labels = labels, lb = lb, ub = ub,
      #   rho1 = rho1, rho2 = rho2
      # )
      fit <- RcppDE::DEoptim(
        lower = lbs,
        upper = ubs,
        fn = EvaluateBessel,
        X = X, labels = labels, lb = lb, ub = ub,
        rho1 = rho1, rho2 = rho2
      )$optim
      fit$par <- fit$bestmem
      # fit <- DEoptimR::JDEoptim(
      #   lower = lbs,
      #   upper = ubs,
      #   fn = EvaluateBessel,
      #   tol = 0.01,
      #   X = X, labels = labels, lb = lb, ub = ub,
      #   rho1 = rho1, rho2 = rho2
      # )
      # fit <- nloptr::stogo(
      #   x0 = x0,
      #   fn = EvaluateBessel,
      #   lower = lbs,
      #   upper = ubs,
      #   X = X, labels = labels, lb = lb, ub = ub,
      #   rho1 = rho1, rho2 = rho2
      # )
    }

    k1 <- fit$par[1]
    alpha1 <- get_alpha(k1, rho1, d)
    k2 <- fit$par[2]
    alpha2 <- get_alpha(k2, rho2, d)
    k12norm <- fit$par[3]
    k12 <- get_k12(k12norm, k1, k2)
    beta12 <- fit$par[4]
    alpha12 <- get_alpha12bis(beta12, alpha1, alpha2)
    tau <- get_tau(k12, alpha12, rho1, rho2, d)

    if (!estimate_rho) {
      return(list(
        rho1 = rho1,
        alpha1 = alpha1,
        rho2 = rho2,
        alpha2 = alpha2,
        tau = tau,
        alpha12 = alpha12,
        fmin = fit$value
      ))
    }
  }

  # Next, estimate also intensities
  # par is now alpha1, alpha2, k12, tau, k1, k2
  if (is.null(init)) {
    x0 <- c(
      alpha1,
      alpha2,
      k12,
      tau,
      get_k(rho1, alpha1, d),
      get_k(rho2, alpha2, d)
    )
  } else {
    x0 <- c(
      init$alpha1,
      init$alpha2,
      get_k(init$tau * sqrt(rho1 * rho2), init$alpha12, d),
      init$tau,
      get_k(rho1, init$alpha1, d),
      get_k(rho2, init$alpha2, d)
    )
  }

  fit <- nloptr::neldermead(
    x0 = x0,
    fn = EvaluateBessel,
    lower = c(
      sqrt(.Machine$double.eps),
      sqrt(.Machine$double.eps),
      0,
      0,
      sqrt(.Machine$double.eps),
      sqrt(.Machine$double.eps)
    ),
    upper = c(
      get_alpha_ub(1 / V, d),
      get_alpha_ub(1 / V, d),
      1,
      1,
      1 - 1e-4,
      1 - 1e-4
    ),
    X = X, labels = labels, lb = lb, ub = ub
  )

  alpha1 <- fit$par[1]
  alpha2 <- fit$par[2]
  k12 <- fit$par[3]
  tau <- fit$par[4]
  k1 <- fit$par[5]
  k2 <- fit$par[6]
  rho1 <- get_rho(k1, alpha1, d)
  rho2 <- get_rho(k2, alpha2, d)
  alpha12 <- get_alpha12(tau, k12, rho1, rho2, d)

  list(
    rho1 = rho1,
    alpha1 = alpha1,
    rho2 = rho2,
    alpha2 = alpha2,
    tau = tau,
    alpha12 = alpha12
  )
}

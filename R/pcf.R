#' Estimation of Stationary Bivariate Bessel DPP via PCF
#'
#' @param X A 2-dimensional marked \code{\link[spatstat]{ppp}} object.
#' @param init A list giving an initial guess at the parameters (default:
#'   \code{NULL}).
#' @param type A character scalar specifying whether PCF contrasts should be
#'   minimized jointly (default) or independently. The latter is faster but
#'   provides biased estimates.
#'
#' @return The estimated model parameters as a list.
#' @export
#'
#' @examples
#' res <- lapply(sim, bessel_pcf_estimation, type = "independent")
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
#'
#' res <- lapply(sim, bessel_pcf_estimation, type = "joint")
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
bessel_pcf_estimation <- function(X, init = NULL, type = "joint") {
  Xs <- spatstat::split.ppp(X)
  d <- length(Xs)
  X1 <- Xs[[1]]
  X2 <- Xs[[2]]
  V <- spatstat::volume(X$window)
  pcf12 <- spatstat::pcfcross(X)
  x12 <- pcf12$iso[22:513]
  r12 <- pcf12$r[22:513]

  # Set rho1 and rho2
  if (is.null(init)) {
    rho1 <- X1$n / V
    rho2 <- X2$n / V
  } else {
    rho1 <- init$rho1
    rho2 <- init$rho2
  }

  if (type == "individual") {
    # First marginal model
    alpha1 <- NULL
    if (!is.null(init)) {
      alpha1 <- ialpha1 <- init$alpha1
      ik1 <- get_k(rho1, ialpha1, d)
    }
    m1 <- fit_marginal_model(X1, V, d, rho1, alpha1)
    rho1 <- m1$rho
    alpha1 <- m1$alpha
    k1 <- m1$k

    # Second marginal model
    alpha2 <- NULL
    if (!is.null(init)) {
      alpha2 <- ialpha2 <- init$alpha2
      ik2 <- get_k(rho2, ialpha2, d)
    }
    m2 <- fit_marginal_model(X2, V, d, rho2, alpha2)
    rho2 <- m2$rho
    alpha2 <- m2$alpha
    k2 <- m2$k

    # Cross contrast for tau and alpha12
    funcontrast <- function(par) {
      k12 <- par[1]
      tau <- par[2]
      alpha12 <- get_alpha12(tau, k12, rho1, rho2, d)
      if (alpha12 <= max(alpha1, alpha2))
        return(1e6)
      sum((x12^0.5 - pcftheocross(c(k12, alpha12), r12, rho1, rho2, d)^0.5)^2)
    }
    if (is.null(init)) {
      lbs <- c(0, 0)
      ubs <- c(
        get_k12_ub(k1, k2),
        1
      )
      fit12 <- nloptr::directL(
        fn = funcontrast,
        lower = lbs,
        upper = ubs,
        original = TRUE
      )
      x0 <- fit12$par
      fit12 <- nloptr::neldermead(
        x0 = x0,
        fn = funcontrast,
        lower = lbs,
        upper = ubs
      )
    } else {
      k12 <- get_k(init$tau * sqrt(rho1 * rho2), init$alpha12, d)
      x0 <- c(k12, init$tau)
      fit12 <- nloptr::neldermead(
        x0 = x0,
        fn = funcontrast,
        lower = c(0, 0),
        upper = c(
          get_k12_ub(ik1, ik2),
          1
        )
      )
    }

    k12 <- fit12$par[1]
    tau <- fit12$par[2]
    alpha12 <- get_alpha12(tau, k12, rho1, rho2, d)

    return(list(
      rho1 = rho1,
      alpha1 = alpha1,
      rho2 = rho2,
      alpha2 = alpha2,
      tau = tau,
      alpha12 = alpha12
    ))
  }

  # Get empirical marginal PCFs
  pcf1 <- spatstat::pcf(X1)
  x1 <- pcf1$iso[22:513]
  r1 <- pcf1$r[22:513]
  pcf2 <- spatstat::pcf(X2)
  x2 <- pcf2$iso[22:513]
  r2 <- pcf2$r[22:513]
  # Get common intersection grid
  rmin <- max(min(r1), min(r2), min(r12))
  rmax <- min(max(r1), max(r2), max(r12))
  rc <- seq(rmin, rmax, length.out = 513-22+1)
  # approx
  x1 <- approx(r1, x1, rc)$y
  x2 <- approx(r2, x2, rc)$y
  x12 <- approx(r12, x12, rc)$y

  # Joint PCF contrast
  joint_pcf <- function(par) {
    alpha1 <- par[1]
    alpha2 <- par[2]
    k12 <- par[3]
    tau <- par[4]
    alpha12 <- get_alpha12(tau, k12, rho1, rho2, d)
    k1 <- get_k(rho1, alpha1, d)
    k2 <- get_k(rho2, alpha2, d)
    if (k12 >= get_k12_ub(k1, k2) | alpha12 < max(alpha1, alpha2))
      return(1e6)
    sum((x1^0.5 - pcftheomarginal(alpha1, rc, d)^0.5)^2) +
          sum((x2^0.5 - pcftheomarginal(alpha2, rc, d)^0.5)^2) +
          sum((x12^0.5 - pcftheocross(c(k12, alpha12), rc, rho1, rho2, d)^0.5)^2)
  }
  lbs <- c(
    sqrt(.Machine$double.eps),
    sqrt(.Machine$double.eps),
    0,
    0
  )
  ubs <- c(
    get_alpha_ub(rho1, d),
    get_alpha_ub(rho2, d),
    1,
    1
  )

  if (is.null(init)) {
    joint_fit <- nloptr::directL(
      fn = joint_pcf,
      lower = lbs,
      upper = ubs,
      original = TRUE
    )
    x0 <- joint_fit$par
  } else {
    x0 <- c(
      init$alpha1,
      init$alpha2,
      get_k(init$tau * sqrt(init$rho1 * init$rho2), init$alpha12, d),
      init$tau
    )
  }

  joint_fit <- nloptr::neldermead(
    x0 = x0,
    fn = joint_pcf,
    lower = lbs,
    upper = ubs
  )

  k12 <- joint_fit$par[3]
  tau <- joint_fit$par[4]
  alpha12 <- get_alpha12(tau, k12, rho1, rho2, d)

  # Output estimated parameters
  list(
    rho1 = rho1,
    alpha1 = joint_fit$par[1],
    rho2 = rho2,
    alpha2 = joint_fit$par[2],
    tau = tau,
    alpha12 = alpha12
  )
}

fit_marginal_model <- function(X, V, d, rho = NULL, alpha = NULL) {
  # Get initial values if not provided
  ## rho
  if (is.null(rho)) rho <- X$n / V
  ## alpha
  alpha_lb <- sqrt(.Machine$double.eps)
  alpha_ub <- get_alpha_ub(rho, d)
  pcf <- spatstat::pcf(X)
  x <- pcf$iso[22:513]
  r <- pcf$r[22:513]
  funcontrast <- function(par) {
    sum((x^0.5 - pcftheomarginal(par, r, d)^0.5)^2)
  }
  if (is.null(alpha)) {
    fit <- nloptr::directL(
      fn = funcontrast,
      lower = alpha_lb,
      upper = alpha_ub,
      original = TRUE
    )
    alpha_init <- fit$par
  } else
    alpha_init <- alpha

  # Final fit
  fit <- nloptr::neldermead(
    x0 = alpha_init,
    fn = funcontrast,
    lower = alpha_lb,
    upper = alpha_ub
  )

  alpha <- fit$par
  k <- get_k(rho, alpha, d)
  list(rho = rho, alpha = alpha, k = k)
}

mfunbessel <- function(r, alpha, d = 2) {
  order <- d / 2
  sapply(r, function(.x) {
    x <- sqrt(2 * d * .x^2) / alpha
    if (x < sqrt(.Machine$double.eps))
      return(1 / gamma(1 + order))
    if (x > 1e5)
      return(2^order * sqrt(2 / pi) * cos(x - alpha * pi / 2 - pi/4) / x^(0.5 + order))
    2^order * besselJ(x, order) / x^order
  })
}

pcftheomarginal <- function(par, r, d = 2, ...) {
  1 - mfunbessel(r, par, d = d)^2
}

pcftheocross <- function(par, r, rho1, rho2, d = 2) {
  k <- par[1]
  alpha <- par[2]
  tau2 <- get_rho(k, alpha, d)^2 / (rho1 * rho2)
  1 - tau2 * mfunbessel(r, alpha, d = d)^2
}

get_k <- function(rho, alpha, d) {
  rho * gamma(d / 2 + 1) * (2 * pi * alpha^2 / d)^(d / 2)
}

get_rho <- function(k, alpha, d) {
  k / (gamma(d / 2 + 1) * (2 * pi * alpha^2 / d)^(d / 2))
}

get_alpha_ub <- function(rho, d = 2) {
  spatstat::dppparbounds(
    model = spatstat::dppBessel(
      lambda = rho,
      sigma = 0,
      d = d
    ),
    name = "alpha"
  )[2]
}

get_k12_ub <- function(k1, k2) {
  sqrt(max(min(k1 * k2, (1 - k1) * (1 - k2)), 0))
}

get_alpha12 <- function(tau, k12, rho1, rho2, d) {
  rho12 <- tau * sqrt(rho1 * rho2)
  if (rho12 < 0.001) return(sqrt(.Machine$double.eps))
  (k12 / (rho12 * (2 * pi / d)^(d / 2) * gamma(1 + d / 2)))^(1 / d)
}

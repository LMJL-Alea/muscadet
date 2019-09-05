#' Estimation of Stationary Bivariate Bessel DPP via PCF
#'
#' @param X A 2-dimensional marked \code{\link[spatstat]{ppp}} object.
#' @param pcftheomarginal A function implementing the theoretical marginal PCF.
#' @param pcftheocross A function implementing the theoretical cross PCF.
#'
#' @return The estimated model parameters as a list.
#' @export
#'
#' @examples
#' tau <- 0.2
#'
#' res <- lapply(sim, bessel_pcf_estimation, tau_init = tau, type = "independent")
#' temp <- unlist(res)
#' par(mfrow = c(2, 2))
#' boxplot(temp[seq(2, 600, 6)], main = 'alpha1')
#' abline(h = 0.03)
#' boxplot(temp[seq(4, 600, 6)], main = 'alpha2')
#' abline(h = 0.03)
#' boxplot(temp[seq(5, 600, 6)], main = 'tau')
#' abline(h = tau)
#' boxplot(temp[seq(6, 600, 6)], main = 'alpha12')
#' abline(h = 0.05)
#'
#' res <- lapply(sim, bessel_pcf_estimation, tau_init = tau, type = "joint")
#' temp <- unlist(res)
#' par(mfrow = c(2, 2))
#' boxplot(temp[seq(2, 600, 6)], main = 'alpha1')
#' abline(h = 0.03)
#' boxplot(temp[seq(4, 600, 6)], main = 'alpha2')
#' abline(h = 0.03)
#' boxplot(temp[seq(5, 600, 6)], main = 'tau')
#' abline(h = tau)
#' boxplot(temp[seq(6, 600, 6)], main = 'alpha12')
#' abline(h = 0.05)
bessel_pcf_estimation <- function(X, tau_init = 0.5, eta = 0.1, type = "joint") {
  Xs <- spatstat::split.ppp(X)
  d <- length(Xs)
  X1 <- Xs[[1]]
  X2 <- Xs[[2]]
  V <- spatstat::volume(X$window)

  # First marginal model
  m1 <- fit_marginal_model(X1, V, d)
  rho1 <- m1$rho
  alpha1 <- m1$alpha
  k1 <- m1$k

  # Second marginal model
  m2 <- fit_marginal_model(X2, V, d)
  rho2 <- m2$rho
  alpha2 <- m2$alpha
  k2 <- m2$k

  # Get a first guess at alpha12 and
  # tau by minimizinf cross-contrast
  ub <- alpha1^d * alpha2^d * min(1, (1 / k1 - 1) * (1 / k2 - 1))
  pcf12 <- spatstat::pcfcross(X)
  x12 <- pcf12$iso[22:513]
  r12 <- pcf12$r[22:513]
  funcontrast <- function(par) {
    tau <- par[1]
    alpha12 <- par[2]
    if (alpha12 < max(alpha1, alpha2) || tau^2 * alpha12^(2 * d) >= ub)
      return(1e6)
    sum((x12^0.5 - pcftheocross(par, r12)^0.5)^2)
  }
  alpha12 <- (1+eta) * max(alpha1, alpha2)
  tau <- tau_init
  while (tau^2 * alpha12^(2*d) >= ub) {
    tau <- tau / sqrt(2)
  }
  fit12 <- stats::optim(
    par = c(tau, alpha12),
    fn = funcontrast,
    lower = c(sqrt(.Machine$double.eps), max(alpha1, alpha2)),
    upper = c(1, Inf),
    method = "L-BFGS-B"
  )
  tau <- fit12$par[1]
  alpha12 <- fit12$par[2]
  independent_fit <- list(
    rho1 = rho1,
    alpha1 = alpha1,
    rho2 = rho2,
    alpha2 = alpha2,
    tau = tau,
    alpha12 = alpha12
  )

  if (type != "joint")
    return(independent_fit)

  # Prepare initialization for joint PCF estimation
  x0 <- c(alpha1, alpha2, alpha12, tau)

  # Get empirical marginal PCFs
  pcf1 <- spatstat::pcf(X1)
  x1 <- pcf1$iso[22:513]
  r1 <- pcf1$r[22:513]
  pcf2 <- spatstat::pcf(X2)
  x2 <- pcf2$iso[22:513]
  r2 <- pcf2$r[22:513]
  # Get common intersection grid
  rmin <- max(r1[1], r2[1], r12[1])
  rmax <- min(r1[length(r1)], r2[length(r2)], r12[length(r12)])
  rc <- seq(rmin, rmax, length.out = 513-22+1)
  # approx
  x1 <- approx(r1, x1, rc)$y
  x2 <- approx(r2, x2, rc)$y
  x12 <- approx(r12, x12, rc)$y

  # Joint PCF contrast
  joint_pcf <- function(par) {
    alpha1 <- par[1]
    alpha2 <- par[2]
    alpha12 <- par[3]
    tau <- par[4]
    k1 <- rho1 * gamma(1 + d/2) * (2 * pi * alpha1^2 / d)^(d/2)
    k2 <- rho2 * gamma(1 + d/2) * (2 * pi * alpha2^2 / d)^(d/2)
    ub <- alpha1^d * alpha2^d * min(1, (1 / k1 - 1) * (1 / k2 - 1))
    if (k1 >= 1 || k2 >= 1 ||
        alpha12 < max(alpha1, alpha2) ||
        tau^2 * alpha12^(2 * d) >= ub)
      return(1e6)
    sum((x1^0.5 - pcftheomarginal(alpha1, rc)^0.5)^2 +
          sum(x2^0.5 - pcftheomarginal(alpha2, rc)^0.5)^2 +
          2 * sum(x12^0.5 - pcftheocross(c(tau, alpha12), rc)^0.5)^2)
  }
  joint_fit <- stats::optim(
    par = x0,
    fn = joint_pcf,
    lower = c(sqrt(.Machine$double.eps), sqrt(.Machine$double.eps), max(alpha1, alpha2), sqrt(.Machine$double.eps)),
    upper = c(sqrt(d / (2*pi)) * (rho1 * gamma(1 + d/2))^(-1/d), sqrt(d / (2*pi)) * (rho2 * gamma(1 + d/2))^(-1/d), Inf, 1),
    method = "L-BFGS-B"
  )

  # Output estimated parameters
  list(
    rho1 = rho1,
    alpha1 = joint_fit$par[1],
    rho2 = rho2,
    alpha2 = joint_fit$par[2],
    tau = joint_fit$par[4],
    alpha12 = joint_fit$par[3]
  )
}

fit_marginal_model <- function(X, V, d) {
  rho <- X$n / V
  alpha_lb <- sqrt(.Machine$double.eps)
  alpha_ub <- spatstat::dppparbounds(
    model = spatstat::dppBessel(lambda = rho, sigma = 0, d = d),
    name = "alpha")[2]
  alpha_init <- (alpha_lb + alpha_ub) / 2
  fit <- spatstat::mincontrast(
    observed = spatstat::pcf(X),
    theoretical = pcftheomarginal,
    startpar = alpha_init,
    ctrl = list(q = 1/2, p = 2, rmin = 0.01),
    method = "Brent",
    lower = alpha_lb,
    upper = alpha_ub
  )
  alpha <- fit$par
  k <- rho * gamma(d / 2 + 1) * (2 * pi * alpha^2 / d)^(d / 2)
  list(rho = rho, alpha = alpha, k = k)
}

mfunbessel <- function(r, alpha, d = 2) {
  order <- d / 2
  x <- sqrt(2 * d * r^2) / alpha
  ifelse(r < sqrt(.Machine$double.eps), 1 / gamma(1 + order), ifelse(x > 1e5, 0, 2^order * besselJ(x, order) / x^order))
}

pcftheomarginal <- function(par, r, d = 2, ...) {
  1 - mfunbessel(r, par)^2
}

pcftheocross <- function(par, r, d = 2) {
  1 - par[1]^2 * mfunbessel(r, par[2])^2
}

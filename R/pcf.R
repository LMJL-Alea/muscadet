contr_marginal <- function(alpha, r, y, p = 0.5, q = 2) {
  yobs <- y^p
  ypred <- (1 - exp(-2 * (r / alpha)^2))^p
  sum(c(0, diff(r)) * abs(yobs - ypred)^q * r)
}

contr_beta <- function(beta, r, y, gamma_max, rho1, rho2) {
  part1 <- gamma_max^3 * beta^3 * exp(-4 * min(r)^2 * beta) / (24 * rho1^2 * rho2^2 * pi^4)
  part2 <- gamma_max^2 * beta^2 / (rho1 * rho2 * pi^2) * sum(c(0, diff(r)) * (1 - y) * exp(-2 * beta * r^2) * r)
  part1 - part2
}

compute_gamma <- function(beta, r, y, gamma_max, rho1, rho2) {
  gamma <- 8 * rho1 * rho2 * pi^2 * exp(4 * beta * min(r)^2) * sum(c(0, diff(r)) * r * (1 - y) * exp(-2 * beta * r^2)) / beta
  if (gamma < 0 | gamma > gamma_max) {
    gamma <- 0
    val1 <- sum(c(0, diff(r)) * (1 - y - gamma * beta^2 * exp(- 2 * beta * r^2) / (rho1 * rho2 * pi^2))^2 * r)
    gamma <- gamma_max
    val2 <- sum(c(0, diff(r)) * (1 - y - gamma * beta^2 * exp(- 2 * beta * r^2) / (rho1 * rho2 * pi^2))^2 * r)
    if (val1 < val2) gamma <- 0
  }
  gamma
}

#' Estimation of Stationary Bivariate 2-dimensional DPP
#'
#' @param X a \code{\link[spatstat]{ppp}} object storing the point pattern.
#' @param rmin_alpha The lower bound on distances that should be taken into
#'   account for estimating marginal alpha parameters (default: index 1).
#' @param rmin_alpha12 The lower bound on distances that should be taken into
#'   account for estimating the crossing alpha parameter (default: index 1).
#' @param rmin_tau The lower bound on distances that should be taken into
#'   account for estimating the correlation (default: index 31).
#' @param tau_min Correlation value used to compute a suitable upper bound for
#'   the crossing alpha parameter (default: 0.1).
#' @param p Power used in the marginal contrasts for estimating \code{alpha1}
#'   and \code{alpha2} (default: 0.2).
#'
#' @return A list with the estimated model parameters in the following order:
#'   \code{rho1}, \code{rho2}, \code{alpha1}, \code{alpha2}, \code{alpha12} and
#'   \code{tau}.
#' @export
#'
#' @examples
#' res <- purrr::map_df(sim_gauss5, estimate)
#' boxplot(res$alpha1)
#' abline(h = 0.03, col = "red")
#' boxplot(res$alpha2)
#' abline(h = 0.03, col = "red")
#' boxplot(res$alpha12)
#' abline(h = 0.035, col = "red")
#' boxplot(res$tau)
#' abline(h = 0.5, col = "red")
estimate <- function(X,
                     model = "Gauss",
                     method = "PCF",
                     rmin_alpha = 1,
                     rmin_alpha12 = 1,
                     rmin_tau = 31,
                     tau_min = 0.1,
                     p = 0.2,
                     divisor_alpha = "d",
                     divisor_alpha12 = "d",
                     bw_alpha = "SJ",
                     bw_alpha12 = "SJ") {
  Xs <- spatstat::split.ppp(X)

  # First estimate marginal intensities
  rho2 <- spatstat::intensity(X)
  rho1 <- as.numeric(rho2[1])
  rho2 <- as.numeric(rho2[2])

  # Estimate alpha1
  pcfemp <- spatstat::pcf(Xs[[1]], bw = bw_alpha, divisor = divisor_alpha)
  alpha_ub <- spatstat::dppparbounds(spatstat::dppGauss(lambda = rho1, d = 2))
  alpha_lb <- alpha_ub[2, 1] + sqrt(.Machine$double.eps)
  alpha_ub <- alpha_ub[2, 2] - sqrt(.Machine$double.eps)
  alpha1 <- optimise(
    f = contr_marginal,
    interval = c(alpha_lb, alpha_ub),
    r = pcfemp$r[rmin_alpha:512],
    y = pcfemp$iso[rmin_alpha:512],
    p = p
  )$minimum

  # Estimate alpha2
  pcfemp <- spatstat::pcf(Xs[[2]], bw = bw_alpha, divisor = divisor_alpha)
  alpha_ub <- spatstat::dppparbounds(spatstat::dppGauss(lambda = rho2, d = 2))
  alpha_lb <- alpha_ub[2, 1] + sqrt(.Machine$double.eps)
  alpha_ub <- alpha_ub[2, 2] - sqrt(.Machine$double.eps)
  alpha2 <- optimise(
    f = contr_marginal,
    interval = c(alpha_lb, alpha_ub),
    r = pcfemp$r[rmin_alpha:512],
    y = pcfemp$iso[rmin_alpha:512],
    p = p
  )$minimum

  # Set bounds for cross parameters
  # We estimate
  # - beta = 1 / alpha12^2
  # - gamma = tau^2 rho1 rho2 pi^2 alpha12^4
  pcfemp <- spatstat::pcfcross(X, bw = bw_alpha12, divisor = divisor_alpha12)
  k1 <- rho1 * pi * alpha1^2 # Model dependent
  k2 <- rho2 * pi * alpha2^2 # Model dependent
  gamma_max <- min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps))
  gamma_max <- max(0, gamma_max)

  # Start by beta estimation
  beta_min <- tau_min / (alpha1 * alpha2)
  beta_max <- 2 / (alpha1^2 + alpha2^2) # Model dependent
  beta <- optimise(
    f = contr_beta,
    interval = c(beta_min, beta_max),
    r = pcfemp$r[rmin_alpha12:512],
    y = pcfemp$iso[rmin_alpha12:512],
    gamma_max = gamma_max,
    rho1 = rho1,
    rho2 = rho2
  )$minimum

  # Retrieve gamma
  gamma <- compute_gamma(
    beta = beta,
    r = pcfemp$r[rmin_tau:512],
    y = pcfemp$iso[rmin_tau:512],
    gamma_max = gamma_max,
    rho1 = rho1,
    rho2 = rho2
  )

  if (method == "PCF")
    return(list(
      rho1 = rho1,
      rho2 = rho2,
      alpha1 = alpha1,
      alpha2 = alpha2,
      alpha12 = 1 / sqrt(beta),
      tau = sqrt(gamma / (rho1 * rho2)) * beta / pi
    ))

  # Code for MLE

}

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
      ifelse(
        alpha12 <= max(alpha1, alpha2),
        1e6,
        sum((x12^0.5 - pcftheocross(c(k12, alpha12), r12, rho1, rho2, d)^0.5)^2)
      )
    }
    if (is.null(init)) {
      lbs <- c(0, 0)
      ubs <- c(
        get_k12_ub(k1, k2),
        1
      )
      # fit12 <- rgenoud::genoud(
      #   fn = funcontrast,
      #   nvars = 2,
      #   Domains = cbind(lbs, ubs),
      #   boundary.enforcement = 1,
      #   print.level = 0
      # )
      # fit12 <- GenSA::GenSA(
      #   fn = funcontrast,
      #   lower = lbs,
      #   upper = ubs
      # )
      # x0 <- fit12$par
      fit12 <- RcppDE::DEoptim(
        fn = funcontrast,
        lower = lbs,
        upper = ubs,
        control = RcppDE::DEoptim.control(trace = FALSE)
      )
      x0 <- fit12$optim$bestmem
      # fit12 <- nloptr::directL(
      #   fn = funcontrast,
      #   lower = lbs,
      #   upper = ubs,
      #   original = TRUE
      # )
      # x0 <- fit12$par
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
      alpha12 = alpha12,
      min_contrast = fit12$value
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
    ifelse(
      k12 >= get_k12_ub(k1, k2) | alpha12 < max(alpha1, alpha2),
      1e6,
      sum((x1^0.5 - pcftheomarginal(alpha1, rc, d)^0.5)^2) +
          sum((x2^0.5 - pcftheomarginal(alpha2, rc, d)^0.5)^2) +
          10 * sum((x12^0.5 - pcftheocross(c(k12, alpha12), rc, rho1, rho2, d)^0.5)^2)
    )
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
    joint_fit <- RcppDE::DEoptim(
      fn = joint_pcf,
      lower = lbs,
      upper = ubs,
      control = RcppDE::DEoptim.control(trace = FALSE)
    )
    x0 <- joint_fit$optim$bestmem
    # joint_fit <- nloptr::directL(
    #   fn = joint_pcf,
    #   lower = lbs,
    #   upper = ubs,
    #   original = TRUE
    # )
    # x0 <- joint_fit$par
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
    alpha_init <- optimise(
      f = funcontrast,
      lower = alpha_lb,
      upper = alpha_ub
    )$minimum
  } else {
    alpha_init <- alpha
  }

  fit <- nloptr::neldermead(
    x0 = alpha_init,
    fn = funcontrast,
    lower = alpha_lb,
    upper = alpha_ub
  )

  alpha <- fit$par
  min_contrast <- fit$value
  k <- get_k(rho, alpha, d)

  list(
    rho = rho,
    alpha = alpha,
    k = k,
    min_contrast = min_contrast
  )
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

pcftheomarginal <- function(alpha, r, d = 2) {
  1 - mfunbessel(r, alpha, d = d)^2
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
  ifelse(rho12 < 0.001, sqrt(.Machine$double.eps), (k12 / (rho12 * (2 * pi / d)^(d / 2) * gamma(1 + d / 2)))^(1 / d))
}

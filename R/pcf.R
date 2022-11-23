#' Estimation of Two-Mark Planar DPPs via the Pair Correlation Function
#'
#' @section Empirical estimation of the pair correlation function:
#' The empirical PCF is computed as a kernel estimate of the PCF in which the
#' contribution from an interpoint distance \eqn{d_{ij}} to the estimate of
#' \eqn{g(r)} is divided:
#' - either by \eqn{r} using optional argument `divisor = "r"` in the functions
#' [spatstat.core::pcf()] and [spatstat.core::pcfcross()];
#' - or by \eqn{d_{ij}} using optional argument `divisor = "d"` in the functions
#' [spatstat.core::pcf()] and [spatstat.core::pcfcross()]; it is intended to
#' improve the bias of the estimator when \eqn{r} is close to zero.
#'
#' @param X An object of class [spatstat.geom::ppp] specifying a planar point
#'   pattern.
#' @param model A string specifying the model to be fitted. Choices are
#'   `"Gauss"` or `"Bessel"`. Defaults to `"Gauss"`.
#' @param method A string specifying the estimation method. Choices are `"PCF"`.
#'   Defaults to `"PCF"`.
#' @param rmin_alpha The lower bound on distances that should be taken into
#'   account for estimating marginal alpha parameters (default: index 2).
#' @param rmin_alpha12 The lower bound on distances that should be taken into
#'   account for estimating the crossing alpha parameter (default: index 2).
#' @param rmin_tau The lower bound on distances that should be taken into
#'   account for estimating the correlation (default: index 2).
#' @param p Power for the distance between empirical and moodel-based PCF
#'   values. Defaults to `2`.
#' @param q Power for pointwise evaluations of the PCF. Defaults to `0.5`.
#' @param divisor_marginal Choice of divisor in the estimation formula. Choices
#'   are `"r"` or `"d"`. See Section Empirical estimation of the pair
#'   correlation function for more details. Defaults to `"d"`.
#' @param divisor_cross Choice of divisor in the estimation formula. Choices are
#'   `"r"` or `"d"`. See Section Empirical estimation of the pair correlation
#'   function for more details. Defaults to `"d"`.
#' @param method A character string specifying the estimation method between
#'   `"profiling"` and `"direct"`. Defaults to `"profiling"`.
#' @param B An integer value specifying the number of samples to be generated in
#'   the bootstrap procedure to approximate the distribution of the test
#'   statistics when testing for absence of correlation between marks. Defaults
#'   to `0L` which does not perform the test at all.
#' @param conf_level A numeric value specifying the confidence level when
#'   testing for absence of correlation between marks. Defaults to `0.95`.
#' @param params A length-4 numeric vector specifying values for the marginal
#'   parameters if known. The order needs to be `rho1`, `rho2`, `alpha1` and
#'   `alpha2`. Defaults to `NULL`, in which case, they are estimated.
#' @param full_bootstrap A boolean specifying whether marginal parameters should
#'   be re-estimated when computing the bootstrapped distribution of the tau
#'   statistic. Defaults to `TRUE`.
#'
#' @return A list with the estimated model parameters in the following order:
#'   `rho1`, `rho2`, `alpha1`, `alpha2`, `k12`, `alpha12` and `tau`. Additional
#'   information pertaining to the test for absence of correlation between marks
#'   are returned in the list as well.
#' @export
#'
#' @examples
#' fit_via_pcf(sim_gauss0[[1]])
fit_via_pcf <- function(X,
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
                        conf_level = 0.95,
                        params = NULL,
                        full_bootstrap = TRUE) {
  divisor_marginal <- match.arg(divisor_marginal, c("d", "r"))
  divisor_cross <- match.arg(divisor_cross, c("d", "r"))

  Xs <- spatstat.geom::split.ppp(X)

  if (is.null(params)) {
    # First estimate marginal intensities
    rho2 <- spatstat.geom::intensity(X)
    rho1 <- as.numeric(rho2[1])
    rho2 <- as.numeric(rho2[2])

    # Estimate alpha1
    alpha1 <- compute_marginal_alpha(
      x = Xs[[1]],
      divisor = divisor_marginal,
      rmin = rmin_alpha,
      q = q,
      p = p,
      model = model
    )

    # Estimate alpha2
    alpha2 <- compute_marginal_alpha(
      x = Xs[[2]],
      divisor = divisor_marginal,
      rmin = rmin_alpha,
      q = q,
      p = p,
      model = model
    )
  } else {
    rho1 <- params[1]
    rho2 <- params[2]
    alpha1 <- params[3]
    alpha2 <- params[4]
  }

  # Get cross PCF for estimating tau and alpha12
  pcfemp <- spatstat.core::pcfcross(X, bw = "SJ", divisor = divisor_cross)

  # Notations: beta = 1 / alpha12^2
  bnds <- get_bounds(
    rho1 = rho1, rho2 = rho2,
    alpha1 = alpha1, alpha2 = alpha2,
    d = 2, model = model
  )
  beta_max <- bnds$beta_max
  k1 <- bnds$k1
  k2 <- bnds$k2
  k12maxSq <- max(0, min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps)))
  # M <- beta_max^2 / (rho1 * rho2 * pi^2) * k12maxSq # TO DO: check
  M <- (alpha1 * alpha2 * beta_max)^2 * min(1, (1 - k1) * (1 - k2) / (k1 * k2)) # tau^2 upper bound

  beta_min <- 0

  if (method == "profiling") {
    .contrast_cross <- function(x) {
      tau2 <- compute_tau2_from_beta(
        beta = x,
        r = pcfemp$r[rmin_tau:length(pcfemp$r)],
        y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
        k1 = k1,
        k2 = k2,
        alpha1 = alpha1,
        alpha2 = alpha2
      )

      # TO DO
      if (tau2 * rho1 * rho2 * pi^2 > k12maxSq * x^2) {
        return(1e6)
      }

      contrast_cross(model = model, beta = x, tau = sqrt(tau2), pcfemp = pcfemp, rmin = rmin_alpha12)
    }

    if (abs(beta_max - beta_min) < .Machine$double.eps^0.5) {
      beta <- beta_max
      fmin <- .contrast_cross(beta_max)
    } else {
      opt <- stats::optimise(
        f = .contrast_cross,
        interval = c(beta_min, beta_max),
        tol = .Machine$double.eps
      )

      opts <- list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
      opt <- nloptr::bobyqa(
        x0 = opt$minimum,
        fn = .contrast_cross,
        lower = beta_min,
        upper = beta_max,
        control = opts
      )
      beta <- opt$par
      fmin <- opt$value
    }

    tau <- sqrt(compute_tau2_from_beta(
      beta = beta,
      r = pcfemp$r[rmin_tau:length(pcfemp$r)],
      y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
      k1 = k1,
      k2 = k2,
      alpha1 = alpha1,
      alpha2 = alpha2
    ))
    k12 <- get_kinf_value(
      rho = tau * sqrt(rho1 * rho2),
      alpha = 1 / sqrt(beta),
      d = 2,
      model = model
    )
  } else {
    if (abs(beta_max - beta_min) < .Machine$double.eps^0.5) {
      .contrast_cross <- function(x) {
        are_params_feasible <- check_parameter_set(
          rho1 = rho1, rho2 = rho2,
          alpha1 = alpha1, alpha2 = alpha2,
          alpha12 = 1 / sqrt(beta_max), tau = x,
          d = 2, model = model
        )
        if (!are_params_feasible)
          return(1e6)
        contrast_cross(model, beta_max, x, pcfemp, rmin_alpha12, q = q, p = p)
      }

      opt <- stats::optimise(
        f = .contrast_cross,
        lower = 0,
        upper = sqrt(M),
        tol = .Machine$double.eps
      )

      opt <- nloptr::bobyqa(
        x0 = opt$minimum,
        fn = .contrast_cross,
        lower = 0,
        upper = sqrt(M),
        control = list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
      )

      beta <- beta_max
      tau <- opt$par
    } else {
      .contrast_cross <- function(x) {
        are_params_feasible <- check_parameter_set(
          rho1 = rho1, rho2 = rho2,
          alpha1 = alpha1, alpha2 = alpha2,
          alpha12 = 1 / sqrt(x[1]), tau = x[2],
          d = 2, model = model
        )
        if (!are_params_feasible)
          return(1e6)
        contrast_cross(model, x[1], x[2], pcfemp, rmin_alpha12, q = q, p = p)
      }
      lower <- c(beta_min, 0)
      upper <- c(beta_max, sqrt(M))

      opt <- nloptr::directL(
        fn = .contrast_cross,
        lower = lower,
        upper = upper,
        control = list(xtol_rel = 1e-6, maxeval = 1e3)
      )

      opt <- nloptr::bobyqa(
        x0 = opt$par,
        fn = .contrast_cross,
        lower = lower,
        upper = upper,
        control = list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
      )

      beta <- opt$par[1]
      tau <- opt$par[2]
    }

    fmin <- opt$value
    if (beta < sqrt(.Machine$double.eps))
      k12 <- sqrt(k12maxSq)
    else
      k12 <- min(get_kinf_value(
        rho = tau * sqrt(rho1 * rho2),
        alpha = 1 / sqrt(beta),
        d = 2, model = model
      ), sqrt(k12maxSq))
  }

  # compute tau first via integrated contrast
  stat_np_obs <- contrast_cross(
    model = model,
    beta = 1,
    tau = 0,
    pcfemp,
    rmin = rmin_tau,
    q = q,
    p = p
  )
  stat_p_obs <- tau
  null_np_distr <- NULL
  null_p_distr <- NULL
  np_reject <- NULL
  p_reject <- NULL
  if (B > 0) {
    cli::cli_alert_info("Testing for absence of correlation...")
    null_values <- compute_bootstrap_stats(
      rho1 = rho1, alpha1 = alpha1,
      rho2 = rho2, alpha2 = alpha2,
      w = X$window,
      B = B,
      model = model,
      rmin_alpha = rmin_alpha,
      rmin_alpha12 = rmin_alpha12,
      rmin_tau = rmin_tau,
      q = q,
      p = p,
      divisor_marginal = divisor_marginal,
      divisor_cross = divisor_cross,
      method = method,
      full_bootstrap = full_bootstrap
    )
    null_np_distr <- null_values$nonparametric
    null_p_distr <- null_values$parametric
    quantile_idx <- round(B * conf_level)
    np_reject <- as.logical(stat_np_obs >= sort(null_np_distr)[quantile_idx])
    p_reject <- as.logical(stat_p_obs >= sort(null_p_distr)[quantile_idx])
  }

  list(
    rho1 = rho1,
    rho2 = rho2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    k12 = k12,
    tau = tau,
    alpha12 = 1 / sqrt(beta),
    fmin = fmin,
    stat_np_obs = stat_np_obs,
    stat_p_obs = stat_p_obs,
    null_np_distr = null_np_distr,
    null_p_distr = null_p_distr,
    np_reject = np_reject,
    p_reject = p_reject
  )
}

compute_bootstrap_stats <- function(rho1, alpha1,
                                    rho2, alpha2,
                                    w = list(xrange = c(0, 1), yrange = c(0, 1)),
                                    B = 100,
                                    model = "Gauss",
                                    rmin_alpha = 2,
                                    rmin_alpha12 = 2,
                                    rmin_tau = 2,
                                    q = 1,
                                    p = 2,
                                    divisor_marginal = "d",
                                    divisor_cross = "d",
                                    method = "profiling",
                                    full_bootstrap = TRUE) {
  w <- spatstat.geom::as.owin(w)
  if (model == "Gauss") {
    m1 <- spatstat.core::dppGauss(lambda = rho1, alpha = alpha1, d = 2)
    m2 <- spatstat.core::dppGauss(lambda = rho2, alpha = alpha2, d = 2)
  } else if (model == "Bessel") {
    m1 <- spatstat.core::dppBessel(lambda = rho1, alpha = alpha1, d = 2)
    m2 <- spatstat.core::dppBessel(lambda = rho2, alpha = alpha2, d = 2)
  } else
    cli::cli_abort("Model {model} is not yet implemented. Currently supported models are {.field Gauss} or {.field Bessel}.")

  data1 <- stats::simulate(m1, nsim = B, w = w)
  data2 <- stats::simulate(m2, nsim = B, w = w)

  boot_data <- purrr::map2(data1, data2, ~ {
    spatstat.geom::ppp(
      x = c(.x$x, .y$x),
      y = c(.x$y, .y$y),
      window = w,
      marks = as.factor(c(rep(1, .x$n), rep(2, .y$n)))
    )
  })

  stat1 <- boot_data |>
    purrr::map(spatstat.core::pcfcross, bw = "SJ", divisor = divisor_cross) |>
    purrr::map_dbl(
      .f = contrast_cross,
      model = model,
      beta = 1, tau = 0,
      rmin = rmin_tau, q = q, p = p
    )

  if (full_bootstrap) {
    stat2 <- boot_data |>
      purrr::map(
        .f = fit_via_pcf,
        model = model,
        rmin_alpha = rmin_alpha,
        rmin_alpha12 = rmin_alpha12,
        rmin_tau = rmin_tau,
        q = q,
        p = p,
        divisor_marginal = divisor_marginal,
        divisor_cross = divisor_cross,
        method = method
      ) |>
      purrr::map_dbl("tau")
  } else {
    stat2 <- boot_data |>
      purrr::map(
        .f = fit_via_pcf,
        model = model,
        rmin_alpha = rmin_alpha,
        rmin_alpha12 = rmin_alpha12,
        rmin_tau = rmin_tau,
        q = q,
        p = p,
        divisor_marginal = divisor_marginal,
        divisor_cross = divisor_cross,
        method = method,
        params = c(rho1, rho2, alpha1, alpha2)
      ) |>
      purrr::map_dbl("tau")
  }

  list(nonparametric = stat1, parametric = stat2)
}

jinc <- function(x, alpha) {
  # J_alpha(x) / (x/2)^alpha * gamma(alpha + 1)
  if (x < sqrt(.Machine$double.eps))
    return(1)
  besselJ(x, alpha) / (x / 2)^alpha * gamma(alpha + 1)
}

contrast_marginal <- function(model, alpha, r, y, q = 0.5, p = 2, d = 2) {
  yobs <- y^q
  ypred <- (1 - get_eta_value(r, 1 / alpha^2, model)^2)^q
  sum(c(0, diff(r)) * abs(yobs - ypred)^p)
}

contrast_cross <- function(model, beta, tau, pcfemp, rmin, q = 1, p = 2) {
  r <- pcfemp$r[rmin:length(pcfemp$r)]
  yobs <- pcfemp$iso[rmin:length(pcfemp$r)]^q
  ypred <- (1 - tau^2 * get_eta_value(r, beta, model)^2)^q
  sum(c(0, diff(r)) * abs(yobs - ypred)^p)
}

compute_tau2_from_beta <- function(beta, r, y, k1, k2, alpha1, alpha2,
                                   model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  eta_val <- get_eta_value(r, beta, model = model)
  I1 <- sum(c(0, diff(r)) * (1 - y) * eta_val^2)
  I2 <- sum(c(0, diff(r)) * eta_val^4)
  tauSq_max <- alpha1 * alpha2 * beta * min(1, (1 - k1) * (1 - k2) / (k1 * k2))
  if (I2 < sqrt(.Machine$double.eps))
    return(tauSq_max)
  tauSq <- I1 / I2
  tauSq <- max(0, tauSq)
  min(tauSq_max, tauSq)
}

compute_marginal_alpha <- function(x, divisor, rmin, q = 0.5, p = 2, model = "Gauss") {
  if (model == "Gauss") {
    alpha_ub <- spatstat.core::dppparbounds(spatstat.core::dppGauss(
      lambda = spatstat.geom::intensity(x),
      d = 2
    ))
  } else if (model == "Bessel") {
    alpha_ub <- spatstat.core::dppparbounds(spatstat.core::dppBessel(
      lambda = spatstat.geom::intensity(x),
      d = 2,
      sigma = 0
    ))
  } else
    cli::cli_abort("Model {model} is not yet implemented. Currently supported models are {.field Gauss} or {.field Bessel}.")

  alpha_lb <- alpha_ub[2, 1] + sqrt(.Machine$double.eps)
  alpha_ub <- alpha_ub[2, 2] - sqrt(.Machine$double.eps)

  pcfemp <- spatstat.core::pcf(x, bw = "SJ", divisor = divisor)

  opt <- stats::optimise(
    f = contrast_marginal,
    interval = c(alpha_lb, alpha_ub),
    model = model,
    r = pcfemp$r[rmin:length(pcfemp$r)],
    y = pcfemp$iso[rmin:length(pcfemp$r)],
    q = q,
    p = p,
    tol = .Machine$double.eps
  )

  opts <- list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
  nloptr::bobyqa(
    x0 = opt$minimum,
    fn = contrast_marginal,
    lower = alpha_lb,
    upper = alpha_ub,
    model = model,
    r = pcfemp$r[rmin:length(pcfemp$r)],
    y = pcfemp$iso[rmin:length(pcfemp$r)],
    q = q,
    p = p,
    control = opts
  )$par
}

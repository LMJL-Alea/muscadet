.validate_parameter_set <- function(rho1, rho2, tau, beta,
                                    alpha1, alpha2, dimension = 2) {
  k11 <- rho1 * alpha1^dimension * pi^(dimension/2)
  k22 <- rho2 * alpha2^dimension * pi^(dimension/2)
  k12sq <- tau^2 * rho1 * rho2 * (pi / beta)^dimension
  k12sq_max <- max(0, min(k11 * k22, (1-k11)*(1-k22)-sqrt(.Machine$double.eps)))
  k11 < 1 && k22 < 1 && k12sq < k12sq_max && beta <= 2 / (alpha1^2 + alpha2^2)
}

contrast_marginal <- function(alpha, r, y, q = 0.5, p = 2) {
  yobs <- y^q
  ypred <- (1 - exp(-2 * (r / alpha)^2))^q
  sum(c(0, diff(r)) * abs(yobs - ypred)^p)
}

contrast_cross <- function(beta, tau, pcfemp, rmin, q = 1, p = 2) {
  r <- pcfemp$r[rmin:length(pcfemp$r)]
  yobs <- pcfemp$iso[rmin:length(pcfemp$r)]^q
  ypred <- (1 - tau^2 * exp(-2 * beta * r^2))^q
  sum(c(0, diff(r)) * abs(yobs - ypred)^p)
}

compute_tau <- function(beta_max, M, r, y) {
  I1 <- sum(c(0, diff(r)) * (1 - y) * (1 - exp(-2*r^2*beta_max)) / (2 * r^2))
  I2 <- sum(c(0, diff(r)) * (1 - exp(-4*r^2*beta_max)) / (4 * r^2))
  tauSq <- I1 / I2
  if (tauSq < 0) return(0)
  if (tauSq > M) return(sqrt(M))
  sqrt(tauSq)
}

compute_tau2_from_beta <- function(beta, r, y, k12SqMax, rho1, rho2) {
  I1 <- sum(c(0, diff(r)) * (1 - y) * exp(-2*beta*r^2))
  I2 <- sum(c(0, diff(r)) * exp(-4*r^2*beta))
  # return(I1/I2)
  tauSq <- I1 / I2
  tauSq <- max(0, tauSq)
  min(k12SqMax * beta^2 / rho1 / rho2 / pi^2, tauSq)
}

combine_bw_marginal <- function(x, divisor, rmin_alpha, q = 0.5, p = 2, bw = c("bcv", "SJ")) {
  alpha_ub <- spatstat.core::dppparbounds(spatstat.core::dppGauss(
    lambda = spatstat.geom::intensity(x),
    d = 2
  ))

  alpha_lb <- alpha_ub[2, 1] + sqrt(.Machine$double.eps)
  alpha_ub <- alpha_ub[2, 2] - sqrt(.Machine$double.eps)

  if (is.null(bw))
    return(.compute_alpha_and_precision(x, bw, divisor, alpha_lb, alpha_ub, p, q, rmin_alpha)$alpha)

  df <- bw %>%
    purrr::map(
      .f = .compute_alpha_and_precision,
      x = x, divisor = divisor,
      alpha_lb = alpha_lb, alpha_ub = alpha_ub,
      p = p, q = q,
      rmin = rmin_alpha
    ) %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  sum(df$precision * df$alpha) / sum(df$precision)
}

.compute_alpha_and_precision <- function(x, bw, divisor, alpha_lb, alpha_ub, p, q, rmin) {
  pcfemp <- spatstat.core::pcf(x, bw = bw, divisor = divisor, var.approx = TRUE)

  list(
    alpha = stats::optimise(
      f = contrast_marginal,
      interval = c(alpha_lb, alpha_ub),
      r = pcfemp$r[rmin:length(pcfemp$r)],
      y = pcfemp$iso[rmin:length(pcfemp$r)],
      q = 0.5,
      p = p,
      tol = .Machine$double.eps
    )$minimum,
    precision = 1 / mean(pcfemp$v[-1])
  )
}

#' Estimation of Stationary Bivariate 2-dimensional DPP
#'
#' @section Empirical estimation of the pair correlation function: The empirical
#'   PCF is computed as a kernel estimate of the PCF in which the contribution
#'   from an interpoint distance $d_{ij}$ to the estimate of $g(r)$ is divided:
#'   - either by $r$ using optional argument `divisor = "r"` in the functions
#'   \code{\link[spatstat.core]{pcf}} and \code{\link[spatstat.core]{pcfcross}};
#'   - or by $d_{ij}$ using optional argument `divisor = "d"` in the functions
#'   \code{\link[spatstat.core]{pcf}} and \code{\link[spatstat.core]{pcfcross}};
#'   it is intended to improve the bias of the estimator when $r$ is close to
#'   zero.
#'
#' @param X A \code{\link[spatstat.geom]{ppp}} object storing the point pattern.
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
#' @param bw_marginal A character vector specifying the bandwidths that should
#'   be used to compute the empirical PCF for estimating the marginal
#'   parameters. Choices are `"NULL"`, `"nrd0"`, `"nrd"`, `"ucv"`, `"bcv"` or
#'   `"SJ"`. Defaults to `NULL` which uses \code{\link[spatstat.core]{pcf}}
#'   default bandwidth.
#' @param bw_cross A string specifying the bandwidth that should be used to
#'   compute the empirical PCF for estimating the cross-type parameters. Choices
#'   are `"NULL"`, `"nrd0"`, `"nrd"`, `"ucv"`, `"bcv"` or `"SJ"`. Defaults to
#'   `NULL` which uses \code{\link[spatstat.core]{pcfcross}} default bandwidth.
#'
#' @return A list with the estimated model parameters in the following order:
#'   `rho1`, `rho2`, `alpha1`, `alpha2`, `k12`, `alpha12` and `tau`.
#' @export
#'
#' @examples
#' #res <- purrr::map_df(sim_gauss5, estimate)
#' #boxplot(res$alpha1)
#' #abline(h = 0.03, col = "red")
#' #boxplot(res$alpha2)
#' #abline(h = 0.03, col = "red")
#' #boxplot(res$alpha12)
#' #abline(h = 0.035, col = "red")
#' #boxplot(res$tau)
#' #abline(h = 0.5, col = "red")
estimate <- function(X,
                     model = "Gauss",
                     rmin_alpha = 2,
                     rmin_alpha12 = 2,
                     rmin_tau = 2,
                     q = 0.5,
                     p = 2,
                     divisor_marginal = "d",
                     divisor_cross = "d",
                     bw_marginal = NULL,
                     bw_cross = NULL,
                     method = "profiling",
                     external_fmin = NULL) {
  divisor_marginal <- match.arg(divisor_marginal, c("d", "r"))
  divisor_cross <- match.arg(divisor_cross, c("d", "r"))
  if (!is.null(bw_marginal))
    bw_marginal <- match.arg(bw_marginal, c("nrd0", "nrd", "ucv", "bcv", "SJ"))
  if (!is.null(bw_cross))
    bw_cross <- match.arg(bw_cross, c("nrd0", "nrd", "ucv", "bcv", "SJ"))

  Xs <- spatstat.geom::split.ppp(X)

  # First estimate marginal intensities
  rho2 <- spatstat.geom::intensity(X)
  rho1 <- as.numeric(rho2[1])
  rho2 <- as.numeric(rho2[2])

  # Estimate alpha1
  # pcfemp <- spatstat::pcf(Xs[[1]], bw = bw_marginal, divisor = divisor_marginal)
  alpha1 <- combine_bw_marginal(
    x = Xs[[1]],
    divisor = divisor_marginal,
    rmin_alpha = rmin_alpha,
    q = q,
    p = p,
    bw = bw_marginal
  )

  # Estimate alpha2
  # pcfemp <- spatstat::pcf(Xs[[2]], bw = bw_marginal, divisor = divisor_marginal)
  alpha2 <- combine_bw_marginal(
    x = Xs[[2]],
    divisor = divisor_marginal,
    rmin_alpha = rmin_alpha,
    q = q,
    p = p,
    bw = bw_marginal
  )

  # Get cross PCF for estimating tau and alpha12
  # pcfemp <- spatstat.core::pcfcross(X, bw = bw_cross, divisor = divisor_cross)
  pcfemp <- spatstat.core::pcfcross(X, bw = "SJ", divisor = "d")
  g <- spatstat.core::pcfcross(X, bw = "SJ", divisor = "r")
  auto_rmin <- which(g$iso - lag(g$iso) > 0)[1]
  slope <- (pcfemp$iso[auto_rmin + 1] - pcfemp$iso[auto_rmin]) / (pcfemp$r[auto_rmin + 1] - pcfemp$r[auto_rmin])
  intercept <- pcfemp$iso[auto_rmin] - slope * pcfemp$r[auto_rmin]
  pcfemp$iso[1:auto_rmin] <- intercept + slope * pcfemp$r[1:auto_rmin]
  # rmin_tau <- which(pcfemp$iso - lag(pcfemp$iso) > 0)[1]

  # Model-dependent variables
  # Notations: beta = 1 / alpha12^2
  beta_max <- 2 / (alpha1^2 + alpha2^2)
  k1 <- rho1 * pi * alpha1^2
  k2 <- rho2 * pi * alpha2^2
  k12maxSq <- max(0, min(k1 * k2, (1-k1)*(1-k2)-sqrt(.Machine$double.eps)))
  M <- beta_max^2 / (rho1 * rho2 * pi^2) * k12maxSq
  # cli::cli_alert_info("Done")

  if (method == "profiling") {
    .contrast_cross <- function(x) {
      tau2 <- compute_tau2_from_beta(
        beta = x,
        r = pcfemp$r[rmin_tau:length(pcfemp$r)],
        y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
        k12SqMax = k12maxSq,
        rho1 = rho1,
        rho2 = rho2
      )
      if (tau2 * rho1 * rho2 * pi^2 / x^2 > k12maxSq)
        return(1e6)
      contrast_cross(beta = x, tau = sqrt(tau2), pcfemp = pcfemp, rmin = rmin_alpha12)
    }

    opt <- stats::optimise(
      f = .contrast_cross,
      interval = c(sqrt(.Machine$double.eps), beta_max),
      tol = .Machine$double.eps
    )
    beta <- opt$minimum
    fmin <- opt$objective

    opts <- list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
    opt <- nloptr::bobyqa(
      x0 = beta,
      fn = .contrast_cross,
      lower = sqrt(.Machine$double.eps),
      upper = beta_max,
      control = opts
    )
    beta <- opt$par
    fmin <- opt$value

    tau <- sqrt(compute_tau2_from_beta(
      beta = beta,
      r = pcfemp$r[rmin_tau:length(pcfemp$r)],
      y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
      k12SqMax = k12maxSq,
      rho1 = rho1,
      rho2 = rho2
    ))
    k12 <- tau * sqrt(rho1 * rho2) * pi / beta
  } else if (method == "marginalization") {
    # compute tau first
    tau <- compute_tau(
      r = pcfemp$r[rmin_tau:length(pcfemp$r)],
      y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
      beta_max = beta_max,
      M = M
    )

    # Notations: k12 = tau sqrt(rho1 rho2) pi alpha12^2
    # = tau * sqrt(rho1 * rho2) * pi / beta
    if (tau < .Machine$double.eps^0.5) {
      k12 <- 0
      fmin <- contrast_cross(1, 0, pcfemp, rmin_alpha12)
    } else if (abs(tau^2 - M) < 1e-4) {
      k12 <- tau * sqrt(rho1 * rho2) * pi / beta_max
      fmin <- contrast_cross(beta_max, sqrt(M), pcfemp, rmin_alpha12)
    } else {
      beta_min <- beta_max * tau / sqrt(M)

      .contrast_cross <- function(beta) {
        contrast_cross(beta, tau, pcfemp, rmin_alpha12)
      }

      opt <- stats::optimise(
        f = .contrast_cross,
        interval = c(beta_min, beta_max),
        tol = .Machine$double.eps
      )
      beta <- opt$minimum
      fmin <- opt$objective

      opts <- list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
      opt <- nloptr::bobyqa(
        x0 = beta,
        fn = .contrast_cross,
        lower = beta_min,
        upper = beta_max,
        control = opts
      )

      beta <- opt$par
      fmin <- opt$value

      k12 <- tau * sqrt(rho1 * rho2) * pi / beta
    }
  } else {
    .contrast_cross <- function(x) {
      if (!.validate_parameter_set(rho1, rho2, x[2], x[1], alpha1, alpha2))
        return(1e6)
      contrast_cross(x[1], x[2], pcfemp, rmin_alpha12)
    }
    lower <- c(0, 0)
    upper <- c(beta_max, sqrt(M))

    if (is.null(external_fmin)) {
      opt <- nloptr::directL(
        fn = .contrast_cross,
        lower = lower,
        upper = upper,
        control = list(xtol_rel = 1e-6, maxeval = 1e3)
      )

      x0 <- opt$par
      opt <- nloptr::bobyqa(
        x0 = x0,
        fn = .contrast_cross,
        lower = lower,
        upper = upper,
        control = list(xtol_rel = .Machine$double.eps, maxeval = 1e6)
      )
    } else {
      x0 <- (lower + upper) / 2
      opt <- nloptr::bobyqa(
        x0 = x0,
        fn = .contrast_cross,
        lower = lower,
        upper = upper,
        control = list(stopval = external_fmin, maxeval = 1e6)
      )
    }

    beta <- opt$par[1]
    tau <- opt$par[2]
    fmin <- opt$value
    if (beta < sqrt(.Machine$double.eps))
      k12 <- sqrt(k12maxSq)
    else
      k12 <- min(tau * sqrt(rho1 * rho2) * pi / beta, sqrt(k12maxSq))
  }

  list(
    rho1 = rho1,
    rho2 = rho2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    k12 = k12,
    tau = tau,
    alpha12 = ifelse(
      tau < 1e-4,
      NA,
      sqrt(k12 / tau / sqrt(rho1 * rho2) / pi)
    ),
    fmin = fmin
  )
}

contrast_marginal <- function(alpha, r, y, q = 0.5, p = 2) {
  yobs <- y^q
  ypred <- (1 - exp(-2 * (r / alpha)^2))^q
  sum(c(0, diff(r)) * abs(yobs - ypred)^p)
}

compute_tau <- function(beta_max, M, r, y) {
  I1 <- sum(c(0, diff(r)) * (1 - y) * (1 - exp(-2*r^2*beta_max)) / (2 * r^2))
  I2 <- sum(c(0, diff(r)) * (1 - exp(-4*r^2*beta_max)) / (4 * r^2))
  tauSq <- I1 / I2
  if (tauSq < 0) return(0)
  if (tauSq > M^2 * beta_max^2) return(M * beta_max)
  sqrt(tauSq)
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
      q = q,
      p = p
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
                     method = "PCF",
                     rmin_alpha = 2,
                     rmin_alpha12 = 2,
                     rmin_tau = 2,
                     q = 0.5,
                     p = 2,
                     divisor_marginal = "d",
                     divisor_cross = "d",
                     bw_marginal = NULL,
                     bw_cross = NULL) {
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
  pcfemp <- spatstat.core::pcfcross(X, bw = bw_cross, divisor = divisor_cross)

  # Model-dependent variables
  # Notations: beta = 1 / alpha12^2
  beta_max <- 2 / (alpha1^2 + alpha2^2)
  k1 <- rho1 * pi * alpha1^2
  k2 <- rho2 * pi * alpha2^2
  k12max <- sqrt(max(0, min(k1 * k2, (1-k1)*(1-k2)-sqrt(.Machine$double.eps))))
  M <- k12max / sqrt(rho1 * rho2) / pi

  # compute tau first
  tau <- compute_tau(
    r = pcfemp$r[rmin_tau:length(pcfemp$r)],
    y = pcfemp$iso[rmin_tau:length(pcfemp$r)],
    beta_max = beta_max,
    M = M
  )

  # Notations: gamma = tau sqrt(rho1 rho2) pi alpha12^2
  gamma_min <- tau * sqrt(rho1 * rho2) * pi / beta_max
  gamma_max <- sqrt(max(0, min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps))))

  if (gamma_max - gamma_min < sqrt(.Machine$double.eps))
    gamma <- gamma_min
  else {
    contrast_cross <- function(gamma) {
      r <- pcfemp$r[rmin_alpha12:length(pcfemp$r)]
      yobs <- pcfemp$iso[rmin_alpha12:length(pcfemp$r)]^q
      ypred <- (1 - tau^2 * exp(-2 * r^2 * tau * sqrt(rho1*rho2) * pi / gamma))^q
      sum(c(0, diff(r)) * abs(yobs - ypred)^p)
    }

    gamma <- stats::optimise(
      f = contrast_cross,
      interval = c(gamma_min, gamma_max)
    )$minimum
  }

  if (method == "PCF")
    return(list(
      rho1 = rho1,
      rho2 = rho2,
      alpha1 = alpha1,
      alpha2 = alpha2,
      k12 = gamma,
      tau = tau,
      alpha12 = ifelse(
        tau < sqrt(.Machine$double.eps),
        NA,
        sqrt(gamma / tau / sqrt(rho1 * rho2) / pi)
      )
    ))

  # Code for MLE

}

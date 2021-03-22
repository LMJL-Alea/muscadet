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
  df <- bw %>%
    purrr::map(~ {
      pcfemp <- spatstat.core::pcf(x, bw = .x, divisor = divisor, var.approx = TRUE)
      c(
        alpha = stats::optimise(
          f = contrast_marginal,
          interval = c(alpha_lb, alpha_ub),
          r = pcfemp$r[rmin_alpha:length(pcfemp$r)],
          y = pcfemp$iso[rmin_alpha:length(pcfemp$r)],
          q = q,
          p = p
        )$minimum,
        precision = 1 / mean(pcfemp$v[-1])
      )
    }) %>%
    purrr::transpose() %>%
    purrr::simplify_all()
  sum(df$precision * df$alpha) / sum(df$precision)
}

#' Estimation of Stationary Bivariate 2-dimensional DPP
#'
#' @param X a \code{\link[spatstat.geom]{ppp}} object storing the point pattern.
#' @param rmin_alpha The lower bound on distances that should be taken into
#'   account for estimating marginal alpha parameters (default: index 1).
#' @param rmin_alpha12 The lower bound on distances that should be taken into
#'   account for estimating the crossing alpha parameter (default: index 1).
#' @param rmin_tau The lower bound on distances that should be taken into
#'   account for estimating the correlation (default: index 31).
#' @param p Power used in the marginal contrasts for estimating \code{alpha1}
#'   and \code{alpha2} (default: 0.2).
#'
#' @return A list with the estimated model parameters in the following order:
#'   \code{rho1}, \code{rho2}, \code{alpha1}, \code{alpha2}, \code{alpha12} and
#'   \code{tau}.
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
    bw = c("bcv", "SJ")
  )

  # Estimate alpha2
  # pcfemp <- spatstat::pcf(Xs[[2]], bw = bw_marginal, divisor = divisor_marginal)
  alpha2 <- combine_bw_marginal(
    x = Xs[[2]],
    divisor = divisor_marginal,
    rmin_alpha = rmin_alpha,
    q = q,
    p = p,
    bw = c("bcv", "SJ")
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

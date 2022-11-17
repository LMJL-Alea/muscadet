get_bounds <- function(rho1, rho2, alpha1, alpha2,
                       d = 2, model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  beta_max <- 1 / get_alpha12_lb(alpha1, alpha2, model)^2
  k1 <- get_kinf_value(rho1, alpha1, d, model)
  k2 <- get_kinf_value(rho2, alpha2, d, model)
  list(beta_max = beta_max, k1 = k1, k2 = k2)
}

get_eta_value <- function(r, beta,
                          d = 2, model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = exp(-beta * r^2),
    Bessel = jinc(sqrt(2 * d * beta) * r, d / 2)
  )
}

get_kinf_value <- function(rho, alpha,
                           d = 2, model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = rho * (pi * alpha^2)^(d / 2),
    Bessel = rho * (2 * pi * alpha^2 / d)^(d / 2) * gamma(d / 2 + 1)
  )
}

get_alpha12_lb <- function(alpha1, alpha2,
                           model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = sqrt((alpha1^2 + alpha2^2) / 2),
    Bessel = max(alpha1, alpha2)
  )
}

get_khat_value <- function(r, rho, alpha,
                           d = 2, model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = {
      alphaSq <- alpha^2
      rho * (pi * alphaSq)^(d / 2) * exp(-pi^2 * alphaSq * r^2)
    },
    Bessel = ifelse(
      2 * pi^2 * r^2 * alpha^2 < d,
      rho * (2 * pi * alpha^2 / d)^(d / 2) * gamma(d / 2 + 1),
      0
    )
  )
}

get_khat_matrix <- function(r, rho1, rho2, alpha1, alpha2, alpha12, tau,
                            d = 2, model = c("Gauss", "Bessel")) {
  model <- rlang::arg_match(model)
  k1 <- get_khat_value(r, rho1, alpha1, d, model)
  k2 <- get_khat_value(r, rho2, alpha2, d, model)
  k12 <- get_khat_value(r, tau * sqrt(rho1 * rho2), alpha12, d, model)
  matrix(c(k1, k12, k12, k2), nrow = 2, ncol = 2)
}

check_parameter_set <- function(rho1, rho2, alpha1, alpha2, alpha12, tau,
                                d = 2, model = c("Gauss", "Bessel")) {
  k1 <- get_kinf_value(rho1, alpha1, d, model)
  if (k1 <= 0 || k1 >= 1)
    return(FALSE)
  k2 <- get_kinf_value(rho2, alpha2, d, model)
  if (k2 <= 0 || k2 >= 1)
    return(FALSE)
  k12 <- get_kinf_value(tau * sqrt(rho1 * rho2), alpha12, d, model)
  if (k12^2 >= min(k1 * k2, (1 - k1) * (1 - k2)))
    return(FALSE)
  if (alpha12 < get_alpha12_lb(alpha1, alpha2, model))
    return(FALSE)
  TRUE
}

#' Feasible range for marginal repulsion rates
#'
#' For a given intensity, dimension and model, provides bounds for the repulsion
#' rate.
#'
#' @param rho A numeric value specifying the intensity of the DPP.
#' @param d An integer value specifying the dimension of the DPP. Defaults to `2L`.
#' @param model A string specifying the model to be used. Choices are `"Gauss"` or `"Bessel"`. Defaults to `"Gauss"`.
#'
#' @return Displays the interval to the console but outputs `NULL`.
#' @export
#'
#' @examples
#' get_feasible_repulsion_rate(rho = 100)
get_feasible_repulsion_rate <- function(rho, d = 2, model = c("Gauss", "Bessel")) {
  alpha_ub <- 1 / sqrt(get_kinf_value(
    rho = rho, alpha = 1,
    d = d, model = model
  ))
  cli::cli_alert_info(" Marginal repulsion rate: (0, {alpha_ub}).")
}

#' Feasible range for cross parameters
#'
#' For a given set of marginal parameters (intensities and repulsion rates),
#' provides bounds for the cross repulsion rate and the correlation between
#' marks.
#'
#' @param rho1 A numeric value specifying the intensity of the first marginal
#'   DPP.
#' @param rho2 A numeric value specifying the intensity of the second marginal
#'   DPP.
#' @param alpha1 A numeric value specifying the repulsion rate of the first
#'   marginal DPP.
#' @param alpha2 A numeric value specifying the repulsion rate of the second
#'   marginal DPP.
#' @param d An integer value specifying the dimension of the DPP. Defaults to
#'   `2L`.
#' @param model A string specifying the model to be used. Choices are `"Gauss"`
#'   or `"Bessel"`. Defaults to `"Gauss"`.
#'
#' @return Displays the intervals to the console but outputs `NULL`.
#' @export
#'
#' @examples
#' get_feasible_space(rho1 = 100, rho2 = 100, alpha1 = 0.03, alpha2 = 0.03)
get_feasible_space <- function(rho1, rho2, alpha1, alpha2,
                               d = 2, model = c("Gauss", "Bessel")) {
  alpha1_ub <- 1 / sqrt(get_kinf_value(rho = rho1, alpha = 1, d = d, model = model))
  if (!is.null(alpha1) && (alpha1 <= 0 || alpha1 >= alpha1_ub)) {
    cli::cli_alert_info(" First marginal repulsion rate: (0, {alpha1_ub}).")
    cli::cli_abort("The input first marginal repulsion rate is not feasible.")
  }
  k1 <- get_kinf_value(rho = rho1, alpha = alpha1, d = d, model = model)

  alpha2_ub <- 1 / sqrt(get_kinf_value(rho = rho2, alpha = 1, d = d, model = model))
  if (!is.null(alpha2) && (alpha2 <= 0 || alpha2 >= alpha2_ub)) {
    cli::cli_alert_info("Second marginal repulsion rate: (0, {alpha2_ub}).")
    cli::cli_abort("The input second marginal repulsion rate is not feasible.")
  }
  k2 <- get_kinf_value(rho = rho2, alpha = alpha2, d = d, model = model)

  alpha12_lb <- get_alpha12_lb(alpha1 = alpha1, alpha2 = alpha2, model = model)
  cli::cli_alert_info("Cross repulsion rate: [{alpha12_lb}, {Inf}).")

  k12 <- get_kinf_value(
    rho = sqrt(rho1 * rho2),
    alpha = alpha12_lb,
    d = d, model = model
  )

  tau_ub <- sqrt(max(0, min(k1 * k2, (1 - k1) * (1 - k2)))) / k12
  cli::cli_alert_info("Correlation: [0, {tau_ub}).")
}

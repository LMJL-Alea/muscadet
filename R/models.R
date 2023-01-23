#' Available Parametric Determinantal Point Processes
#'
#' Lists available DPP parametric models in the package.
#'
#' @return A character vector storing the available parametric DPP models.
#' @export
#'
#' @examples
#' MUSCADET_DPP_MODELS()
MUSCADET_DPP_MODELS <- function() {
  c("Gauss", "Bessel")
}

jinc <- function(x, alpha) {
  # J_alpha(x) / (x/2)^alpha * gamma(alpha + 1)
  if (x < sqrt(.Machine$double.eps))
    return(1)
  besselJ(x, alpha) / (x / 2)^alpha * gamma(alpha + 1)
}

get_bounds <- function(rho1, rho2, alpha1, alpha2,
                       d = 2, model = MUSCADET_DPP_MODELS()) {
  model <- rlang::arg_match(model)
  beta_max <- 1 / get_alpha12_lb(alpha1, alpha2, model)^2
  k1 <- get_kinf_value(rho1, alpha1, d, model)
  k2 <- get_kinf_value(rho2, alpha2, d, model)
  list(beta_max = beta_max, k1 = k1, k2 = k2)
}

get_eta_value <- function(r, beta,
                          d = 2, model = MUSCADET_DPP_MODELS()) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = exp(-beta * r^2),
    Bessel = jinc(sqrt(2 * d * beta) * r, d / 2)
  )
}

get_kinf_value <- function(rho, alpha,
                           d = 2, model = MUSCADET_DPP_MODELS()) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = rho * (pi * alpha^2)^(d / 2),
    Bessel = rho * (2 * pi * alpha^2 / d)^(d / 2) * gamma(d / 2 + 1)
  )
}

get_alpha12_lb <- function(alpha1, alpha2,
                           model = MUSCADET_DPP_MODELS()) {
  model <- rlang::arg_match(model)
  switch(
    model,
    Gauss = sqrt((alpha1^2 + alpha2^2) / 2),
    Bessel = max(alpha1, alpha2)
  )
}

get_khat_value <- function(r, rho, alpha,
                           d = 2, model = MUSCADET_DPP_MODELS()) {
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
                            d = 2, model = MUSCADET_DPP_MODELS()) {
  model <- rlang::arg_match(model)
  k1 <- get_khat_value(r, rho1, alpha1, d, model)
  k2 <- get_khat_value(r, rho2, alpha2, d, model)
  k12 <- get_khat_value(r, tau * sqrt(rho1 * rho2), alpha12, d, model)
  matrix(c(k1, k12, k12, k2), nrow = 2, ncol = 2)
}

check_parameter_set <- function(rho1, rho2, alpha1, alpha2, alpha12, tau,
                                d = 2, model = MUSCADET_DPP_MODELS()) {
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

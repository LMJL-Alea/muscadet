#' Estimation of multi-mark DPP via Maximum Likelihood
#'
#' @param data An object of class [spatstat.geom::ppp] specifying an observed
#'   planar determinantal point process.
#' @param initial_guess A numeric vector specifying an initial guess for the
#'   model parameters that maximize the likelihood. Defaults to `NULL`, which
#'   initializes at \eqn{0.5} all parameters after suitable transformation into
#'   \eqn{[0, 1]}. If provided, expected order is `c(alpha1, alpha2, alpha12,
#'   tau)`.
#' @param fixed_marginal_parameters A boolean value specifying whether the
#'   marginal parameters should be estimated separately using the marginal
#'   likelihood and then fixed to further estimate the cross parameters.
#'   Defaults to `FALSE`.
#' @param model A string specifying a DPP model. Choices are `"gauss"` or
#'   `"bessel"`. Defaults to `"gauss"`.
#' @param optimizer A string specifying one of the available derivative-free
#'   optimizers in NLOpt. Defaults to `"bobyqa"`.
#' @param num_threads An integer value specifying the number of thread to run
#'   on. Defaults to `1L`.
#' @param N An integer value specifying the maximum truncation index for Fourier
#'   transform. Defaults to `512L`.
#' @param verbose_level An integer value specifying the display information
#'   during optimization. Choose `0L` for no information, `1L` for global
#'   information at likelihood setup or `2L` for detailed information at each
#'   function evaluation. Defaults to `0L.
#'
#' @return A list with two components:
#' - `par` (optimized model parameters);
#' - `value` (`-2 logLik` at the maximum likelihood).
#'
#' @export
#' @examples
#' opt <- mle(sim_gauss5[[1]])
mle <- function(data,
                initial_guess = NULL,
                fixed_marginal_parameters = FALSE,
                model = "gauss",
                optimizer = "bobyqa",
                num_threads = 1,
                N = 512,
                verbose_level = 0) {
  points <- cbind(data$x, data$y)
  marks <- data$marks

  # Bivariate case
  # Fixed marginal params so need to estimate them through marginal likelihood
  marginal_parameters <- NULL
  if (!is.null(marks) && fixed_marginal_parameters) {
    # Initial point available, use it
    if (!is.null(initial_guess)) {
      alpha1 <- initial_guess[["alpha1"]]
      alpha2 <- initial_guess[["alpha2"]]
      init <- c(alpha1, alpha2)
      marginal_parameters <- data |>
        spatstat.geom::split.ppp() |>
        purrr::map2(
          .y = init,
          .f = mle,
          model = model,
          optimizer = optimizer,
          num_threads = num_threads,
          N = N,
          verbose_level = verbose_level
        ) |>
        purrr::map("par") |>
        purrr::map_dbl("alpha")
    } else { # Naive initialization
      marginal_parameters <- data |>
        spatstat.geom::split.ppp() |>
        purrr::map(
          .f = mle,
          model = model,
          optimizer = optimizer,
          num_threads = num_threads,
          N = N,
          verbose_level = verbose_level
        ) |>
        purrr::map("par") |>
        purrr::map_dbl("alpha")
    }
  }

  lower_bound <- c(data$window$xrange[1], data$window$yrange[1])
  upper_bound <- c(data$window$xrange[2], data$window$yrange[2])

  nd_grid <- generate_nd_grid(N, dim(points)[2], lower_bound, upper_bound)

  if (!is.null(marks) && !is.null(initial_guess)) {
    if (fixed_marginal_parameters) {
      rho1 <- initial_guess[["rho1"]]
      rho2 <- initial_guess[["rho2"]]
      alpha1 <- marginal_parameters[1]
      alpha2 <- marginal_parameters[2]
      k12 <- initial_guess[["k12"]]
      tau <- initial_guess[["tau"]]
      initial_guess <- NULL
      if (validate_init(rho1, rho2, alpha1, alpha2, k12, tau, 2))
        initial_guess <- c(k12, tau)
      else
        message("Unfeasible initial guess")
    }
    else {
      alpha1 <- initial_guess[["alpha1"]]
      alpha2 <- initial_guess[["alpha2"]]
      k12 <- initial_guess[["k12"]]
      tau <- initial_guess[["tau"]]
      initial_guess <- c(alpha1, alpha2, k12, tau)
    }
  }

  dpp <- methods::new(DeterminantalPointProcess)
  dpp$SetLikelihoodModel(model);
  dpp$SetOptimizer(optimizer);
  dpp$Fit(
    points, lower_bound, upper_bound, nd_grid,
    initial_guess, marks, marginal_parameters,
    num_threads, N, verbose_level
  )
}

validate_init <- function(rho1, rho2, alpha1, alpha2, k12, tau, dimension) {
  k1 <- rho1 * alpha1^dimension * pi^(dimension / 2)
  k2 <- rho2 * alpha2^dimension * pi^(dimension / 2)
  alpha12_min <- sqrt((alpha1 * alpha1 + alpha2 * alpha2) / 2) # Gaussian-specific
  k12_min <- tau * sqrt(rho1 * rho2) * alpha12_min^dimension * pi^(dimension / 2) # Gaussian-specific
  k12_ub <- sqrt(max(0, min(k1 * k2, (1-k1)*(1-k2)-sqrt(.Machine$double.eps))))
  if (k12 > k12_ub) {
    message("k12 exceeds upper bound")
    print(k12)
    print(k12_ub)
    return(FALSE)
  }
  if (k12_min > k12_ub) {
    message("alpha12 is below lower bound")
    return(FALSE)
  }
  TRUE
}

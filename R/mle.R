#' Estimation of multi-mark DPP via Maximum Likelihood
#'
#' @param data A \code{\link[spatstat]{ppp}} object storing an observed planar
#'   determinantal point process.
#' @param initial_guess A vector providing an intial guess for the model
#'   parameters that maximize the likelihood. Defaults to \code{NULL}, which
#'   initializes at \code{0.5} all parameters after suitable transformation into
#'   [0,1]. If provided, expected order is `c(alpha1, alpha2, alpha12, tau)`.
#' @param fixed_marginal_parameters Boolean specifying whether the marginal
#'   parameters should be estimated separately using the marginal likelihood and
#'   then fixed to further estimate the cross parameters. Default is `FALSE`.
#' @param model A DPP model. For now either gauss or bessel (default: gauss).
#' @param optimizer Any derivative-free optimizer in NlOpt (default: bobyqa).
#' @param num_threads Number of threas to run on (default: 1).
#' @param N Maximum truncation index for Fourier transform (default: 512).
#' @param verbose_level Display information during optimization. Choose \code{0}
#'   for no information, \code{1} for global information at likelihood setup or
#'   \code{2} for detailed  information at each function evaluation. Defaults to
#'   \code{0}.
#'
#' @return A list with two components: \code{par} (optimized model parameters)
#'   and \code{value} (\code{-2 logLik} at the maximum likelihood).
#' @export
#'
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
      marginal_parameters <- data %>%
        spatstat.geom::split.ppp() %>%
        purrr::map2(
          .y = init,
          .f = mle,
          model = model,
          optimizer = optimizer,
          num_threads = num_threads,
          N = N,
          verbose_level = verbose_level
        ) %>%
        purrr::map("par") %>%
        purrr::map_dbl("alpha")
    } else { # Naive initialization
      marginal_parameters <- data %>%
        spatstat.geom::split.ppp() %>%
        purrr::map(
          .f = mle,
          model = model,
          optimizer = optimizer,
          num_threads = num_threads,
          N = N,
          verbose_level = verbose_level
        ) %>%
        purrr::map("par") %>%
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

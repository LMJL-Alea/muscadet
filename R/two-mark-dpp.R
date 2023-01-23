#' Two-Mark Determinantal Point Process Model
#'
#' @description Implements a generic two-mark planar determinantal point process
#'   model. It is possible to set its parameters, fit the model to an observed
#'   planar point pattern or simulate a planar point pattern.
#'
#' @param value A numeric value specifying the value of a given parameter of the
#'   model.
#'
#' @export
TwoMarkDPP <- R6::R6Class(
  classname = "TwoMarkDPP",
  public = list(
    #' @description The [TwoMarkDPP] constructor.
    #'
    #' @param model A string specifying the underlying parametric model to be
    #'   used. Choices are `"Gauss"` or `"Bessel"`. Defaults to `"Gauss"`.
    #'
    #' @return The [TwoMarkDPP] object invisibly.
    #'
    #' @examples
    #' mod <- TwoMarkDPP$new(model = "Gauss")
    initialize = function(model = MUSCADET_DPP_MODELS) {
      model <- rlang::arg_match(model)
      private$m_ParametricModel <- model
    },

    #' @description Estimates the parameters of the two-mark DPP from an
    #'   observed planar point pattern.
    #'
    #' @param X An object of class [spatstat.geom::ppp].
    #' @param method A string specifying the kind of estimator to use. Choices
    #'   are either `"pcf"` or `"mle"`. Defaults to `"pcf"`.
    #'
    #' @return The [TwoMarkDPP] object invisibly.
    #'
    #' @examples
    #' mod <- TwoMarkDPP$new(model = "Gauss")
    #' mod$fit(sim_gauss0[[1]])
    fit = function(X, method = c("pcf", "mle")) {
      if (!inherits(X, "ppp"))
        cli::cli_abort("The input {.arg X} should be of class {.cls ppp}.")
      method <- rlang::arg_match(method)
      private$m_WindowSize <- X$window$xrange[2]
      if (X$window$yrange[2] != private$m_WindowSize)
        cli::cli_abort("Currently only supports square windows.")
      params <- switch(
        method,
        pcf = fit_via_pcf(X, model = private$m_ParametricModel),
        mle = fit_via_mle(X, model = private$m_ParametricModel)
      )
      private$m_FirstIntensity <- params[["rho1"]]
      private$m_SecondIntensity <- params[["rho2"]]
      private$m_FirstRepulsionRate <- params[["alpha1"]]
      private$m_SecondRepulsionRate <- params[["alpha2"]]
      private$m_CrossRepulsionRate <- params[["alpha12"]]
      private$m_BetweenMarkCorrelation <- params[["tau"]]

      private$m_UpperFirstIntensity <- 1 / get_kinf_value(
        rho = 1,
        alpha = private$m_FirstRepulsionRate,
        model = private$m_ParametricModel
      )
      private$m_UpperSecondIntensity <- 1 / get_kinf_value(
        rho = 1,
        alpha = private$m_SecondRepulsionRate,
        model = private$m_ParametricModel
      )
      private$m_UpperFirstRepulsionRate <- sqrt(1 / get_kinf_value(
        rho = private$m_FirstIntensity,
        alpha = 1,
        model = private$m_ParametricModel
      ))
      private$m_UpperSecondRepulsionRate <- sqrt(1 / get_kinf_value(
        rho = private$m_SecondIntensity,
        alpha = 1,
        model = private$m_ParametricModel
      ))
      private$m_LowerCrossRepulsionRate <- get_alpha12_lb(
        alpha1 = private$m_FirstRepulsionRate,
        alpha2 = private$m_SecondRepulsionRate,
        model = private$m_ParametricModel
      )
      k11 <- get_kinf_value(
        rho = private$m_FirstIntensity,
        alpha = private$m_FirstRepulsionRate,
        model = private$m_ParametricModel
      )
      k22 <- get_kinf_value(
        rho = private$m_SecondIntensity,
        alpha = private$m_SecondRepulsionRate,
        model = private$m_ParametricModel
      )
      k12 <- get_kinf_value(
        rho = sqrt(private$m_FirstIntensity * private$m_SecondIntensity),
        alpha = private$m_CrossRepulsionRate,
        model = private$m_ParametricModel
      )
      private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12

      invisible(self)
    },

    #' @description Samples one or more planar point patterns according to the
    #'   two-mark DPP.
    #'
    #' @param n An integer value specifying the sample size.
    #'
    #' @return A [base::list] of objects of class [spatstat.geom::ppp] storing
    #'   sampled planar point patterns.
    #'
    #' @examples
    #' mod <- TwoMarkDPP$new(model = "Gauss")
    #' mod$fit(sim_gauss0[[1]])
    #' mod$random(n = 1)
    random = function(n) {
      if (is.null(private$m_FirstIntensity) ||
          is.null(private$m_SecondIntensity) ||
          is.null(private$m_FirstRepulsionRate) ||
          is.null(private$m_SecondRepulsionRate) ||
          is.null(private$m_CrossRepulsionRate) ||
          is.null(private$m_BetweenMarkCorrelation) ||
          is.null(private$m_WindowSize))
        cli::cli_abort("Some parameters of the two-mark DPP have not been set or
                       estimated.")
      rbidpp(
        n = n,
        rho1 = private$m_FirstIntensity,
        rho2 = private$m_SecondIntensity,
        alpha1 = private$m_FirstRepulsionRate,
        alpha2 = private$m_SecondRepulsionRate,
        alpha12 = private$m_CrossRepulsionRate,
        tau = private$m_BetweenMarkCorrelation,
        L = private$m_WindowSize,
        model = private$m_ParametricModel
      )
    },

    #' @description Prints the [TwoMarkDPP] object.
    #'
    #' @param ... Extra parameters to be passed on to next methods. Current
    #'   unused.
    #'
    #' @return NULL
    #'
    #' @examples
    #' mod <- TwoMarkDPP$new()
    #' mod
    print = function(...) {
      rho1 <- private$m_FirstIntensity
      if (is.null(rho1)) {
        if (is.null(private$m_UpperFirstIntensity))
          rho1 <- "?"
        else
          rho1 <- paste0("(", 0, ", ", round(private$m_UpperFirstIntensity, digits = 3), ")")
      } else rho1 <- round(rho1, digits = 3)

      rho2 <- private$m_SecondIntensity
      if (is.null(rho2)) {
        if (is.null(private$m_UpperSecondIntensity))
          rho2 <- "?"
        else
          rho2 <- paste0("(", 0, ", ", round(private$m_UpperSecondIntensity, digits = 3), ")")
      } else rho2 <- round(rho2, digits = 3)

      alpha1 <- private$m_FirstRepulsionRate
      if (is.null(alpha1)) {
        if (is.null(private$m_UpperFirstRepulsionRate))
          alpha1 <- "?"
        else
          alpha1 <- paste0("(", 0, ", ", round(private$m_UpperFirstRepulsionRate, digits = 3), ")")
      } else alpha1 <- round(alpha1, digits = 3)

      alpha2 <- private$m_SecondRepulsionRate
      if (is.null(alpha2)) {
        if (is.null(private$m_UpperSecondRepulsionRate))
          alpha2 <- "?"
        else
          alpha2 <- paste0("(", 0, ", ", round(private$m_UpperSecondRepulsionRate, digits = 3), ")")
      } else alpha2 <- round(alpha2, digits = 3)

      alpha12 <- private$m_CrossRepulsionRate
      if (is.null(alpha12)) {
        if (is.null(private$m_LowerCrossRepulsionRate))
          alpha12 <- "?"
        else
          alpha12 <- paste0("[", round(private$m_LowerCrossRepulsionRate, digits = 3), ", ", Inf, ")")
      } else alpha12 <- round(alpha12, digits = 3)

      tau <- private$m_BetweenMarkCorrelation
      if (is.null(tau)) {
        if (is.null(private$m_UpperBetweenMarkCorrelation))
          tau <- "?"
        else
          tau <- paste0("[", 0, ", ", round(private$m_UpperBetweenMarkCorrelation, digits = 3), ")")
      } else tau <- round(tau, digits = 3)

      L <- private$m_WindowSize
      if (is.null(L)) L <- "?"

      cli::cli_h1("Two-mark planar {private$m_ParametricModel} determinantal point process")
      cli::cli_bullets(c(
        "*" = "Intensity of the 1st marginal DPP: {rho1}",
        "*" = "Intensity of the 2nd marginal DPP: {rho2}",
        "*" = "Repulsion rate of the 1st marginal DPP: {alpha1}",
        "*" = "Repulsion rate of the 2nd marginal DPP: {alpha2}",
        "*" = "Repulsion rate between marks: {alpha12}",
        "*" = "Correlation between marks: {tau}",
        "*" = "Window size: {L}"
      ))
    }
  ),
  active = list(
    #' @field first_intensity Setter and getter for the first intensity
    #'   \eqn{\rho_1}.
    first_intensity = function(value) {
      if (missing(value)) return(private$m_FirstIntensity)
      if (value <= 0)
        cli::cli_abort("The intensity of the 1st marginal DPP cannot be
                         negative.")
      if (!is.null(private$m_UpperFirstIntensity)) {
        if (value >= private$m_UpperFirstIntensity)
          cli::cli_abort("The intensity of the 1st marginal DPP cannot be
                           greater than {round(private$m_UpperFirstIntensity,
                           digits = 3)}.")
      }
      private$m_FirstIntensity <- value
      private$m_UpperFirstRepulsionRate <- sqrt(1 / get_kinf_value(
        rho = value,
        alpha = 1,
        model = private$m_ParametricModel
      ))
      if (!is.null(private$m_SecondIntensity) &&
          !is.null(private$m_FirstRepulsionRate) &&
          !is.null(private$m_SecondRepulsionRate) &&
          !is.null(private$m_CrossRepulsionRate)) {
        k11 <- get_kinf_value(
          rho = value,
          alpha = private$m_FirstRepulsionRate,
          model = private$m_ParametricModel
        )
        k22 <- get_kinf_value(
          rho = private$m_SecondIntensity,
          alpha = private$m_SecondRepulsionRate,
          model = private$m_ParametricModel
        )
        k12 <- get_kinf_value(
          rho = sqrt(value * private$m_SecondIntensity),
          alpha = private$m_CrossRepulsionRate,
          model = private$m_ParametricModel
        )
        private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12
      }
    },

    #' @field second_intensity Setter and getter for the second intensity
    #'   \eqn{\rho_2}.
    second_intensity = function(value) {
      if (missing(value)) return(private$m_SecondIntensity)
      if (value <= 0)
        cli::cli_abort("The intensity of the 2nd marginal DPP cannot be
                         negative.")
      if (!is.null(private$m_UpperSecondIntensity)) {
        if (value >= private$m_UpperSecondIntensity)
          cli::cli_abort("The intensity of the 2nd marginal DPP cannot be
                           greater than {round(private$m_UpperSecondIntensity,
                           digits = 3)}.")
      }
      private$m_SecondIntensity <- value
      private$m_UpperSecondRepulsionRate <- sqrt(1 / get_kinf_value(
        rho = value,
        alpha = 1,
        model = private$m_ParametricModel
      ))
      if (!is.null(private$m_FirstIntensity) &&
          !is.null(private$m_FirstRepulsionRate) &&
          !is.null(private$m_SecondRepulsionRate) &&
          !is.null(private$m_CrossRepulsionRate)) {
        k11 <- get_kinf_value(
          rho = private$m_FirstIntensity,
          alpha = private$m_FirstRepulsionRate,
          model = private$m_ParametricModel
        )
        k22 <- get_kinf_value(
          rho = value,
          alpha = private$m_SecondRepulsionRate,
          model = private$m_ParametricModel
        )
        k12 <- get_kinf_value(
          rho = sqrt(private$m_FirstIntensity * value),
          alpha = private$m_CrossRepulsionRate,
          model = private$m_ParametricModel
        )
        private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12
      }
    },

    #' @field first_repulsion_rate Setter and getter for the first repulsion
    #'   rate \eqn{\alpha_1}.
    first_repulsion_rate = function(value) {
      if (missing(value)) return(private$m_FirstRepulsionRate)
      if (value <= 0)
        cli::cli_abort("The repulsion rate of the 1st marginal DPP cannot be
                         negative.")
      if (!is.null(private$m_UpperFirstRepulsionRate)) {
        if (value >= private$m_UpperFirstRepulsionRate)
          cli::cli_abort("The repulsion rate of the 1st marginal DPP cannot be
                           greater than {round(private$m_UpperFirstRepulsionRate,
                           digits = 3)}.")
      }
      private$m_FirstRepulsionRate <- value
      private$m_UpperFirstIntensity <- 1 / get_kinf_value(
        rho = 1,
        alpha = value,
        model = private$m_ParametricModel
      )
      if (!is.null(private$m_SecondRepulsionRate)) {
        private$m_LowerCrossRepulsionRate <- get_alpha12_lb(
          alpha1 = value,
          alpha2 = private$m_SecondRepulsionRate,
          model = private$m_ParametricModel
        )
        if (!is.null(private$m_FirstIntensity) &&
            !is.null(private$m_SecondIntensity) &&
            !is.null(private$m_CrossRepulsionRate)) {
          k11 <- get_kinf_value(
            rho = private$m_FirstIntensity,
            alpha = value,
            model = private$m_ParametricModel
          )
          k22 <- get_kinf_value(
            rho = private$m_SecondIntensity,
            alpha = private$m_SecondRepulsionRate,
            model = private$m_ParametricModel
          )
          k12 <- get_kinf_value(
            rho = sqrt(private$m_FirstIntensity * private$m_SecondIntensity),
            alpha = private$m_CrossRepulsionRate,
            model = private$m_ParametricModel
          )
          private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12
        }
      }
    },

    #' @field second_repulsion_rate Setter and getter for the second repulsion
    #'   rate \eqn{\alpha_2}.
    second_repulsion_rate = function(value) {
      if (missing(value)) return(private$m_SecondRepulsionRate)
      if (value <= 0)
        cli::cli_abort("The repulsion rate of the 2nd marginal DPP cannot be
                         negative.")
      if (!is.null(private$m_UpperSecondRepulsionRate)) {
        if (value >= private$m_UpperSecondRepulsionRate)
          cli::cli_abort("The repulsion rate of the 2nd marginal DPP cannot be
                           greater than {round(private$m_UpperSecondRepulsionRate,
                           digits = 3)}.")
      }
      private$m_SecondRepulsionRate <- value
      private$m_UpperSecondIntensity <- 1 / get_kinf_value(
        rho = 1,
        alpha = value,
        model = private$m_ParametricModel
      )
      if (!is.null(private$m_FirstRepulsionRate)) {
        private$m_LowerCrossRepulsionRate <- get_alpha12_lb(
          alpha1 = private$m_FirstRepulsionRate,
          alpha2 = value,
          model = private$m_ParametricModel
        )
        if (!is.null(private$m_FirstIntensity) &&
            !is.null(private$m_SecondIntensity) &&
            !is.null(private$m_CrossRepulsionRate)) {
          k11 <- get_kinf_value(
            rho = private$m_FirstIntensity,
            alpha = private$m_FirstRepulsionRate,
            model = private$m_ParametricModel
          )
          k22 <- get_kinf_value(
            rho = private$m_SecondIntensity,
            alpha = value,
            model = private$m_ParametricModel
          )
          k12 <- get_kinf_value(
            rho = sqrt(private$m_FirstIntensity * private$m_SecondIntensity),
            alpha = private$m_CrossRepulsionRate,
            model = private$m_ParametricModel
          )
          private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12
        }
      }
    },

    #' @field cross_repulsion_rate Setter and getter for the cross repulsion
    #'   rate \eqn{\alpha_{12}}.
    cross_repulsion_rate = function(value) {
      if (missing(value)) return(private$m_CrossRepulsionRate)
      if (!is.null(private$m_LowerCrossRepulsionRate)) {
        if (value < private$m_LowerCrossRepulsionRate)
          cli::cli_abort("The repulsion rate between marks cannot be lower than
                         {round(private$m_LowerCrossRepulsionRate, digits = 3)}.
                         ")
      }
      private$m_CrossRepulsionRate <- value
      if (!is.null(private$m_FirstIntensity) &&
          !is.null(private$m_SecondIntensity) &&
          !is.null(private$m_FirstRepulsionRate) &&
          !is.null(private$m_SecondRepulsionRate)) {
        k11 <- get_kinf_value(
          rho = private$m_FirstIntensity,
          alpha = private$m_FirstRepulsionRate,
          model = private$m_ParametricModel
        )
        k22 <- get_kinf_value(
          rho = private$m_SecondIntensity,
          alpha = private$m_SecondRepulsionRate,
          model = private$m_ParametricModel
        )
        k12 <- get_kinf_value(
          rho = sqrt(private$m_FirstIntensity * private$m_SecondIntensity),
          alpha = value,
          model = private$m_ParametricModel
        )
        private$m_UpperBetweenMarkCorrelation <- sqrt(max(0, min(k11 * k22, (1 - k11) * (1 - k22)))) / k12
      }
    },

    #' @field between_mark_correlation Setter and getter for the correlation
    #'   \eqn{\tau} between marks.
    between_mark_correlation = function(value) {
      if (missing(value)) return(private$m_BetweenMarkCorrelation)
      if (value < 0)
        cli::cli_abort("The correlation between marks cannot be negative.")
      if (!is.null(private$m_UpperBetweenMarkCorrelation)) {
        if (value >= private$m_UpperBetweenMarkCorrelation)
          cli::cli_abort("The correlation between marks cannot be greater than
                         {round(private$m_UpperBetweenMarkCorrelation, digits =
                         3)}.")
      }
      private$m_BetweenMarkCorrelation <- value
    },

    #' @field window_size Setter and getter for the window size \eqn{L}.
    window_size = function(value) {
      if (missing(value)) return(private$m_WindowSize)
      private$m_WindowSize <- value
    }
  ),
  private = list(
    m_FirstIntensity = NULL,
    m_UpperFirstIntensity = NULL,

    m_SecondIntensity = NULL,
    m_UpperSecondIntensity = NULL,

    m_FirstRepulsionRate = NULL,
    m_UpperFirstRepulsionRate = NULL,

    m_SecondRepulsionRate = NULL,
    m_UpperSecondRepulsionRate = NULL,

    m_CrossRepulsionRate = NULL,
    m_LowerCrossRepulsionRate = NULL,

    m_BetweenMarkCorrelation = NULL,
    m_UpperBetweenMarkCorrelation = NULL,

    m_WindowSize = NULL,
    m_ParametricModel = NULL
  )
)

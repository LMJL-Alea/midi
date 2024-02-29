#' Cylinder radial compartment class
#'
#' @description A class to model restricted diffusion in a cylinder in the plane
#'   perpendicular to the cylinder axis.
#'
#' @param radius A numeric value specifying the radius of the cylinder in
#'   meters.
#' @param radius_sd A numeric value specifying the standard deviation of the
#'   radius of the cylinder in meters. Defaults to `NULL`. If specified, a Gamma
#'   distribution of axon radii is used with shape and scale parameters computed
#'   so that the mean is equal to `radius` and the variance is equal to
#'   `radius_sd^2`.
CylinderRadialCompartment <- R6::R6Class(
  "CylinderRadialCompartment",
  public = list(
    #' @description Instantiates a new cylinder radial compartment.
    #'
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   cylinder in m\eqn{^2}.s\eqn{^{-1}}.
    #'
    #' @return An instance of the [`CylinderRadialCompartment`] class.
    initialize = function(radius, diffusivity, radius_sd = NULL) {
      self$set_radius(radius, radius_sd)
      private$diffusivity <- diffusivity
    },

    #' @description Sets the radius of the cylinder.
    #'
    #' @return The instance of the [`CylinderRadialCompartment`] class,
    #'   invisibly.
    set_radius = function(radius, radius_sd = NULL) {
      private$radius <- radius
      private$integrate_over_radius <- FALSE
      if (!is.null(radius_sd)) {
        private$radius_sd <- radius_sd
        private$integrate_over_radius <- TRUE
        radius_var <- radius_sd^2
        shape <- radius^2 / radius_var
        scale <- radius_var / radius
        withr::with_seed(1234, {
          private$radius_sample <- stats::rgamma(
            n = 1000L,
            shape = shape,
            scale = scale
          )
        })
      }
      invisible(self)
    },

    #' @description Computes the signal attenuation predicted by the model.
    #'
    #' @param small_delta A numeric value specifying the duration of the
    #'  gradient pulse in seconds.
    #' @param big_delta A numeric value specifying the duration between the
    #'  gradient pulses in seconds.
    #' @param G A numeric value specifying the strength of the gradient in
    #'  T.m\eqn{^{-1}}.
    #' @param echo_time A numeric value specifying the echo time in seconds.
    #' @param n_max An integer value specifying the maximum order of the Bessel
    #'  function. Defaults to `20L`.
    #' @param m_max An integer value specifying the maximum number of extrema
    #' for the Bessel function. Defaults to `50L`.
    #'
    #' @return A numeric value storing the predicted signal attenuation.
    #'
    #' @examples
    #' sodermanComp <- SodermanCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9
    #' )
    #' sodermanComp$get_signal(0.03, 0.03, 0.040)
    #'
    #' staniszComp <- StaniszCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9
    #' )
    #' staniszComp$get_signal(0.03, 0.03, 0.040)
    #'
    #' neumanComp <- NeumanCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9
    #' )
    #' neumanComp$get_signal(0.03, 0.03, 0.040, echo_time = 0.040)
    #'
    #' callaghanComp <- CallaghanCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9
    #' )
    #' callaghanComp$get_signal(0.03, 0.03, 0.040)
    #'
    #' vanGelderenComp <- VanGelderenCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9
    #' )
    #' vanGelderenComp$get_signal(0.03, 0.03, 0.040)
    get_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      if (private$integrate_over_radius) {
        private$compute_integrated_signal(
          small_delta = small_delta,
          big_delta = big_delta,
          G = G,
          echo_time = echo_time,
          n_max = n_max,
          m_max = m_max
        )
      } else {
        private$compute_signal(
          small_delta = small_delta,
          big_delta = big_delta,
          G = G,
          echo_time = echo_time,
          n_max = n_max,
          m_max = m_max
        )
      }
    }
  ),
  private = list(
    radius = NULL,
    diffusivity = NULL,
    radius_sd = NULL,
    radius_sample = NULL,
    integrate_over_radius = FALSE,
    gamma = 2.675e8, # rad s^-1 T^-1,
    besselJ_derivative = function(x, n) {
      ifelse(
        n == 0,
        -besselJ(x, 1),
        (besselJ(x, n - 1) - besselJ(x, n + 1)) / 2
      )
    },
    compute_integrated_signal = function(small_delta,
                                         big_delta,
                                         G,
                                         echo_time,
                                         n_max,
                                         m_max) {
      work_compartment <- self$clone()
      fun <- function(x) {
        work_compartment$set_radius(x)
        work_compartment$get_signal(
          small_delta = small_delta, big_delta = big_delta, G = G,
          echo_time = echo_time, n_max = n_max, m_max = m_max
        )
      }
      mean(sapply(private$radius_sample, fun))
    }
  )
)

#' Soderman's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Soderman's model.
#'
#' @export
SodermanCompartment <- R6::R6Class(
  "SodermanCompartment",
  inherit = CylinderRadialCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      x <- private$gamma * small_delta * G * private$radius
      if (x < .Machine$double.eps) return(1)
      (2 * besselJ(x, 1) / x)^2
    }
  )
)

#' Stanisz's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Stanisz's model.
#'
#' @export
StaniszCompartment <- R6::R6Class(
  "StaniszCompartment",
  inherit = CylinderRadialCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      x <- private$gamma * small_delta * G * private$radius
      if (x < .Machine$double.eps) return(1)
      n <- seq_len(n_max)
      work_vector <- exp(-n^2 * pi^2 * private$diffusivity * big_delta / private$radius^2)
      work_vector <- work_vector * (1 - (-1)^n * cos(x))
      work_vector <- work_vector / ((x^2 - (n * pi)^2)^2)
      2 * (1 - cos(x)) / x^2 + 4 * x^2 * sum(work_vector)
    }
  )
)

#' Neuman's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Neuman's model.
#'
#' @export
NeumanCompartment <- R6::R6Class(
  "NeumanCompartment",
  inherit = CylinderRadialCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      if (is.null(echo_time)) {
        cli::cli_abort("The echo time must be specified
                       through the {.arg echo_time} argument.")
      }
      first_term <- 7 * private$gamma^2 * small_delta^2 * G^2 * private$radius^4 / (48 * private$diffusivity * echo_time)
      second_term <- 2 - 99 * private$radius^2 / (56 * private$diffusivity * echo_time)
      second_term <- max(0, second_term)
      exp(-first_term * second_term)
    }
  )
)

#' Callaghan's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Callaghan's model.
#'
#' @export
CallaghanCompartment <- R6::R6Class(
  "CallaghanCompartment",
  inherit = CylinderRadialCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      orders <- 0:n_max
      extrms <- bessel_extrema[orders + 1, 1:m_max]
      work_value <- private$gamma * small_delta * G * private$radius
      bessel_value <- private$besselJ_derivative(work_value, orders)
      first_term <- 4 * 2^(orders != 0)
      second_term <- exp(-extrms^2 * private$diffusivity * big_delta / private$radius^2)
      third_term <- 1 / (1 - (orders / extrms)^2)
      third_term[1, 1] <- 1
      fourth_term <- (work_value * bessel_value / (work_value^2 - extrms^2))^2
      if (work_value < sqrt(.Machine$double.eps)) fourth_term[1, 1] <- 1 / 4
      sum(first_term * second_term * third_term * fourth_term)
    }
  )
)

#' Van Gelderen's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Van Gelderen's model.
#'
#' @export
VanGelderenCompartment <- R6::R6Class(
  "VanGelderenCompartment",
  inherit = CylinderRadialCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      alpha_m_sq <- (bessel_extrema[2, 1:m_max] / private$radius)^2
      K <- -2 * private$gamma^2 * G^2
      first_term <- 2 * private$diffusivity * alpha_m_sq * small_delta - 2
      second_term <- 2 * (exp(-private$diffusivity * alpha_m_sq * small_delta) + exp(-private$diffusivity * alpha_m_sq * big_delta))
      third_term <- exp(-private$diffusivity * alpha_m_sq * (big_delta - small_delta)) + exp(-private$diffusivity * alpha_m_sq * (big_delta + small_delta))
      denom_term <- private$diffusivity^2 * alpha_m_sq^3 * (private$radius^2 * alpha_m_sq - 1)
      sum_terms <- (first_term + second_term - third_term) / denom_term
      exp(sum(sum_terms) * K)
    }
  )
)

#' Circularly-shaped compartment class
#'
#' @description A class to model restricted diffusion in a bounded medium
#'   described by a circular shape with a given radius.
#'
#' @keywords internal
CircularlyShapedCompartment <- R6::R6Class(
  "CircularlyShapedCompartment",
  inherit = BaseCompartment,
  public = list(
    #' @description Instantiates a new circular compartment.
    #'
    #' @param radius A numeric value specifying the radius of the sphere in
    #'   \eqn{\mu}m.
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   sphere in \eqn{\mu}m\eqn{^2}.ms\eqn{^{-1}}.
    #'
    #' @return An instance of the [`CircularlyShapedCompartment`] class.
    initialize = function(radius = NULL, diffusivity = NULL) {
      # Set the radius
      if (is.null(radius))
        private$radius <- default_sphere_radius()
      else {
        if (radius <= 0)
          cli::cli_abort("The sphere radius should be strictly positive.")
        private$radius <- radius
      }

      # Set the diffusivity
      if (is.null(diffusivity))
        private$diffusivity <- default_free_diffusivity()
      else {
        if (diffusivity <= 0)
          cli::cli_abort("The sphere diffusivity should be strictly positive.")
        private$diffusivity <- diffusivity
      }
    }
  ),
  private = list(
    radius = NULL,
    diffusivity = NULL,
    parameters = function() list(
      Radius = private$radius,
      Diffusivity = private$diffusivity
    )
  )
)

#' Soderman's model for restricted diffusion in a cylinder
#'
#' @description A class to model restricted diffusion in a cylinder using the
#'   Soderman's model.
#'
#' @references Söderman, O., & Jönsson, B. (1995). Restricted diffusion in
#'   cylindrical geometry. Journal of Magnetic Resonance, Series A, 117(1),
#'   94-97.
#'
#' @export
SodermanCompartment <- R6::R6Class(
  "SodermanCompartment",
  inherit = CircularlyShapedCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      x <- gyromagnetic_ratio() * small_delta * G * private$radius
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
#' @references Stanisz, G. J., Wright, G. A., Henkelman, R. M., & Szafer, A.
#'   (1997). An analytical model of restricted diffusion in bovine optic nerve.
#'   Magnetic Resonance in Medicine, 37(1), 103-111.
#'
#' @export
StaniszCompartment <- R6::R6Class(
  "StaniszCompartment",
  inherit = CircularlyShapedCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      x <- gyromagnetic_ratio() * small_delta * G * private$radius
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
#' @references Neuman, C. H. (1974). Spin echo of spins diffusing in a bounded
#'   medium. The Journal of Chemical Physics, 60(11), 4508-4511.
#'
#' @export
NeumanCompartment <- R6::R6Class(
  "NeumanCompartment",
  inherit = CircularlyShapedCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      if (is.null(echo_time)) {
        cli::cli_abort("The echo time must be specified
                       through the {.arg echo_time} argument.")
      }
      first_term <- 7 * gyromagnetic_ratio()^2 * small_delta^2 * G^2 * private$radius^4 / (48 * private$diffusivity * echo_time)
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
#' @references Callaghan, P. T. (1995). Pulsed-gradient spin-echo NMR for
#'   planar, cylindrical, and spherical pores under conditions of wall
#'   relaxation. Journal of magnetic resonance, Series A, 113(1), 53-59.
#'
#' @export
CallaghanCompartment <- R6::R6Class(
  "CallaghanCompartment",
  inherit = CircularlyShapedCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      orders <- 0:n_max
      extrms <- bessel_extrema[orders + 1, 1:m_max]
      work_value <- gyromagnetic_ratio() * small_delta * G * private$radius
      bessel_value <- besselJ_derivative(work_value, orders)
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
#' @references Vangelderen, P., DesPres, D., Vanzijl, P. C. M., & Moonen, C. T.
#'   W. (1994). Evaluation of restricted diffusion in cylinders. Phosphocreatine
#'   in rabbit leg muscle. Journal of Magnetic Resonance, Series B, 103(3),
#'   255-260.
#'
#' @export
VanGelderenCompartment <- R6::R6Class(
  "VanGelderenCompartment",
  inherit = CircularlyShapedCompartment,
  private = list(
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      alpha_m_sq <- (bessel_extrema[2, 1:m_max] / private$radius)^2
      K <- -2 * gyromagnetic_ratio()^2 * G^2
      first_term <- 2 * private$diffusivity * alpha_m_sq * small_delta - 2
      second_term <- 2 * (exp(-private$diffusivity * alpha_m_sq * small_delta) + exp(-private$diffusivity * alpha_m_sq * big_delta))
      third_term <- exp(-private$diffusivity * alpha_m_sq * (big_delta - small_delta)) + exp(-private$diffusivity * alpha_m_sq * (big_delta + small_delta))
      denom_term <- private$diffusivity^2 * alpha_m_sq^3 * (private$radius^2 * alpha_m_sq - 1)
      sum_terms <- (first_term + second_term - third_term) / denom_term
      exp(sum(sum_terms) * K)
    }
  )
)

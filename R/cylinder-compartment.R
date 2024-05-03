#' Cylinder compartment class
#'
#' @description A class to model restricted diffusion in a cylinder.
#'
#' @export
#' @examples
#' cylComp <- CylinderCompartment$new()
#' cylComp$get_signal(small_delta = 30, big_delta = 30, G = 0.040)
#' cylComp$get_parameter_names()
#' cylComp$get_parameters()
CylinderCompartment <- R6::R6Class(
  "CylinderCompartment",
  inherit = CircularlyShapedCompartment,
  public = list(
    #' @description Instantiates a new cylinder compartment.
    #'
    #' @param axis A length-3 numeric vector specifying the axis of the
    #'  cylinder.
    #' @param restricted_compartment An instance of the
    #'   [`CircularlyShapedCompartment`] class specifying the restricted
    #'   compartment within the sphere. Defaults to a Van Gelderen compartment.
    #'
    #' @return An instance of the [`CylinderCompartment`] class.
    initialize = function(axis = c(0, 0, 1),
                          restricted_compartment = VanGelderenCompartment$new()) {
      if (!is.numeric(axis) || length(axis) != 3) {
        cli::cli_abort("The axis must be a length-3 numeric vector.")
      }
      if (abs(sum(axis^2) - 1) > .Machine$double.eps) {
        cli::cli_abort("The axis must be a unit vector.")
      }
      private$axis <- axis

      if (!inherits(restricted_compartment, "CircularlyShapedCompartment")) {
        cli::cli_abort("The restricted compartment must be one of the circularly
                       shaped compartments.")
      }

      params <- restricted_compartment$get_parameters()
      super$initialize(
        radius = params$Radius,
        diffusivity = params$Diffusivity
      )
      private$internal_compartment <- restricted_compartment
    }
  ),
  private = list(
    axis = NULL,
    internal_compartment = NULL,
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      cos_theta_sq <- sum(direction * private$axis)^2
      radial_signal <- private$internal_compartment$get_signal(
        small_delta = small_delta,
        big_delta = big_delta,
        G = G * sqrt(1 - cos_theta_sq),
        direction = direction,
        echo_time = echo_time,
        n_max = n_max,
        m_max = m_max
      )
      b <- bvalue(small_delta, big_delta, G)
      axial_signal <- exp(-b * cos_theta_sq * private$diffusivity)
      radial_signal * axial_signal
    },
    parameters = function() list(
      CylinderAxis = private$axis,
      CylinderRadius = private$radius,
      CylinderDiffusivity = private$diffusivity
    )
  )
)

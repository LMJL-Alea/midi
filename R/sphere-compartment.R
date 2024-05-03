#' Sphere compartment class
#'
#' @description A class to model restricted diffusion in a sphere.
#'
#' @export
#' @examples
#' sphComp <- SphereCompartment$new()
#' sphComp$get_signal(small_delta = 30, big_delta = 30, G = 0.040)
#' sphComp$get_parameter_names()
#' sphComp$get_parameters()
SphereCompartment <- R6::R6Class(
  "SphereCompartment",
  inherit = CircularlyShapedCompartment,
  public = list(
    #' @description Instantiates a new sphere compartment.
    #'
    #' @param restricted_compartment An instance of the
    #'   [`CircularlyShapedCompartment`] class specifying the restricted
    #'   compartment within the sphere. Defaults to a Van Gelderen compartment.
    #'
    #' @return An instance of the [`SphereCompartment`] class.
    initialize = function(restricted_compartment = VanGelderenCompartment$new()) {
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
    internal_compartment = NULL,
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      private$internal_compartment$get_signal(
        small_delta = small_delta,
        big_delta = big_delta,
        G = G,
        direction = direction,
        echo_time = echo_time,
        n_max = n_max,
        m_max = m_max
      )
    },
    parameters = function() list(
      SphereRadius = private$radius,
      SphereDiffusivity = private$diffusivity
    )
  )
)

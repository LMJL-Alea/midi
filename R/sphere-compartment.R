#' Sphere compartment class
#'
#' @description A class to model restricted diffusion in a sphere.
#'
#' @export
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
      param_names <- restricted_compartment$get_parameter_names()
      param_values <- restricted_compartment$get_parameter_values()
      super$initialize(
        radius = param_values[param_names == "Radius"],
        diffusivity = param_values[param_names == "Diffusivity"]
      )
      private$internal_compartment <- restricted_compartment
    }
  ),
  private = list(
    internal_compartment = NULL,
    compute_signal = function(small_delta, big_delta, G,
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      private$internal_compartment$get_signal(
        small_delta, big_delta, G, echo_time, n_max, m_max
      )
    },
    parameter_names = function() c("SphereRadius", "SphereDiffusivity")
  )
)

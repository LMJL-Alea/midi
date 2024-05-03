#' Free compartment class
#'
#' @description A class to model free unconstrained diffusion.
#'
#' @export
FreeCompartment <- R6::R6Class(
  "FreeCompartment",
  inherit = BaseCompartment,
  public = list(
    #' @description Instantiates a new free compartment.
    #'
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   sphere in \eqn{\mu}m\eqn{^2}.ms\eqn{^{-1}}. Defaults to `NULL`, in
    #'   which case the default free diffusivity of 3
    #'   \eqn{\mu}m\eqn{^2}.ms\eqn{^{-1}} is used.
    #'
    #' @return An instance of the [`FreeCompartment`] class.
    initialize = function(diffusivity = NULL) {
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
    diffusivity = NULL,
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
                              echo_time = NULL,
                              n_max = 20L,
                              m_max = 50L) {
      b <- bvalue(small_delta, big_delta, G)
      exp(-b * private$diffusivity)
    },
    parameters = function() list(
      FreeDiffusivity = private$diffusivity
    )
  )
)

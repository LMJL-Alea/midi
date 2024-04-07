#' Sphere compartment class
#'
#' @description A class to model restricted diffusion in a sphere.
#'
#' @export
SphereCompartment <- R6::R6Class(
  "SphereCompartment",
  public = list(
    #' @description Instantiates a new sphere compartment.
    #'
    #' @param radius A numeric value specifying the radius of the sphere in
    #'   meters.
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   sphere in m\eqn{^2}.s\eqn{^{-1}}.
    #' @param model A character string specifying the radial model to use.
    #'   Choices are `"soderman"`, `"callaghan"`, `"stanisz"`, `"neuman"`, and
    #'   `"vangelderen"`. Defaults to `"soderman"`.
    #'
    #' @return An instance of the [`SphereCompartment`] class.
    initialize = function(radius, diffusivity,
                          model = c("soderman", "callaghan", "stanisz",
                                    "neuman", "vangelderen")) {
      model <- rlang::arg_match(model)
      private$internal_compartment <- switch(
        model,
        soderman = SodermanCompartment$new(radius, diffusivity),
        callaghan = CallaghanCompartment$new(radius, diffusivity),
        stanisz = StaniszCompartment$new(radius, diffusivity),
        neuman = NeumanCompartment$new(radius, diffusivity),
        vangelderen = VanGelderenCompartment$new(radius, diffusivity)
      )
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
    #' sphereComp <- SphereCompartment$new(
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9,
    #'   model = "soderman"
    #' )
    #' sphereComp$get_signal(
    #'   small_delta = 0.03,
    #'   big_delta = 0.03,
    #'   G = 0.040
    #' )
    get_signal = function(small_delta, big_delta, G,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      private$internal_compartment$get_signal(
        small_delta = small_delta,
        big_delta = big_delta,
        G = G,
        echo_time = echo_time,
        n_max = n_max,
        m_max = m_max
      )
    }
  ),
  private = list(
    internal_compartment = NULL
  )
)

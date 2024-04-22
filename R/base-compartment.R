#' Base compartment class
#'
#' @description The base class for compartment models.
BaseCompartment <- R6::R6Class(
  "BaseCompartment",
  public = list(
    #' @description Computes the signal attenuation predicted by the model.
    #'
    #' @param small_delta A numeric value specifying the duration of the
    #'  gradient pulse in milliseconds.
    #' @param big_delta A numeric value specifying the duration between the
    #'  gradient pulses in milliseconds.
    #' @param G A numeric value specifying the strength of the gradient in
    #'  mT.\eqn{\mu}m\eqn{^{-1}}.
    #' @param echo_time A numeric value specifying the echo time in
    #'   milliseconds.
    #' @param n_max An integer value specifying the maximum order of the Bessel
    #'  function. Defaults to `20L`.
    #' @param m_max An integer value specifying the maximum number of extrema
    #' for the Bessel function. Defaults to `50L`.
    #'
    #' @return A numeric value storing the predicted signal attenuation.
    #'
    #' @examples
    #' freeComp <- FreeCompartment$new()
    #' freeComp$get_signal(small_delta = 30, big_delta = 30, G = 0.040e-3)
    #'
    #' sphereComp <- SphereCompartment$new()
    #' sphereComp$get_signal(small_delta = 30, big_delta = 30, G = 0.040e-3)
    #'
    #' sodermanComp <- SodermanCompartment$new()
    #' sodermanComp$get_signal(small_delta = 30, big_delta = 30, G = 0.040e-3)
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
      private$compute_signal(
        small_delta = small_delta,
        big_delta = big_delta,
        G = G,
        echo_time = echo_time,
        n_max = n_max,
        m_max = m_max
      )
    },

    #' @description Returns the names of the compartment parameters
    #'
    #' @return A character vector storing the names of the compartment
    #'   parameters.
    #'
    #' @examples
    #' freeComp <- FreeCompartment$new()
    #' freeComp$get_parameter_names()
    get_parameter_names = function() {
      private$parameter_names()
    },

    #' @description Returns the values of the compartment parameters
    #'
    #' @return A numeric vector storing the values of the compartment
    #'   parameters.
    #'
    #' @examples
    #' freeComp <- FreeCompartment$new()
    #' freeComp$get_parameter_values()
    get_parameter_values = function() {
      private$parameter_values()
    }
  )
)

#' Cylinder compartment class
#'
#' @description A class to model restricted diffusion in a cylinder.
#'
#' @export
CylinderCompartment <- R6::R6Class(
  "CylinderCompartment",
  public = list(
    #' @description Instantiates a new cylinder compartment.
    #'
    #' @param axis A length-3 numeric vector specifying the axis of the
    #'  cylinder.
    #' @param radius A numeric value specifying the radius of the cylinder in
    #'   meters.
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   cylinder in m\eqn{^2}.s\eqn{^{-1}}.
    #' @param radial_model A character string specifying the radial model to
    #'   use. Choices are `"soderman"`, `"callaghan"`, `"stanisz"`, `"neuman"`,
    #'   and `"vangelderen"`. Defaults to `"soderman"`.
    #'
    #' @return An instance of the [`CylinderCompartment`] class.
    initialize = function(axis, radius, diffusivity,
                          radial_model = c("soderman", "callaghan", "stanisz",
                                           "neuman", "vangelderen")) {
      radial_model <- rlang::arg_match(radial_model)
      if (abs(sum(axis^2) - 1) > .Machine$double.eps) {
        cli::cli_abort("The axis must be a unit vector.")
      }
      private$axis <- axis
      if (diffusivity <= .Machine$double.eps) {
        cli::cli_abort("The diffusivity must be a positive value.")
      }
      private$diffusivity <- diffusivity
      private$radial_compartment <- switch(
        radial_model,
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
    #' @param direction A length-3 numeric vector specifying the direction of
    #'   the gradient.
    #' @param echo_time A numeric value specifying the echo time in seconds.
    #' @param n_max An integer value specifying the maximum order of the Bessel
    #'  function. Defaults to `20L`.
    #' @param m_max An integer value specifying the maximum number of extrema
    #' for the Bessel function. Defaults to `50L`.
    #'
    #' @return A numeric value storing the predicted signal attenuation.
    #'
    #' @examples
    #' cylinderComp <- CylinderCompartment$new(
    #'   axis = c(0, 0, 1),
    #'   radius = 1e-6,
    #'   diffusivity = 2.0e-9,
    #'   radial_model = "soderman"
    #' )
    #' cylinderComp$get_signal(
    #'   small_delta = 0.03,
    #'   big_delta = 0.03,
    #'   G = 0.040,
    #'   direction = c(0, 0, 1)
    #' )
    get_signal = function(small_delta, big_delta, G, direction,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      cos_theta_sq <- sum(direction * private$axis)^2
      radial_signal <- private$radial_compartment$get_signal(
        small_delta = small_delta,
        big_delta = big_delta,
        G = G * sqrt(1 - cos_theta_sq),
        echo_time = echo_time,
        n_max = n_max,
        m_max = m_max
      )
      b <- bvalue(small_delta, big_delta, G)
      axial_signal <- exp(-b * cos_theta_sq * private$diffusivity)
      radial_signal * axial_signal
    }
  ),
  private = list(
    axis = NULL,
    diffusivity = NULL,
    radial_compartment = NULL
  )
)

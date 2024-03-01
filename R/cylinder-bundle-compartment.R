#' Cylinder bundle compartment class
#'
#' @description A class to model restricted diffusion in a bundle of cylinders.
#'
#' @export
CylinderBundleCompartment <- R6::R6Class(
  "CylinderBundleCompartment",
  public = list(
    #' @description Instantiates a new cylinder bundle compartment.
    #'
    #' @param axis A length-3 numeric vector specifying the mean axis of the
    #'  cylinder population.
    #' @param radius A numeric value specifying the mean radius of the cylinder
    #'   population in meters.
    #' @param diffusivity A numeric value specifying the diffusivity within the
    #'   cylinders in m\eqn{^2}.s\eqn{^{-1}}.
    #' @param axial_diffusivity A numeric value specifying the axial diffusivity
    #'   in the space outside the cylinders in m\eqn{^2}.s\eqn{^{-1}}.
    #' @param radial_diffusivity A numeric value specifying the radial
    #'   diffusivity in the space outside the cylinders in
    #'   m\eqn{^2}.s\eqn{^{-1}}.
    #' @param cylinder_density A numeric value specifying the density of the
    #'  cylinders in the voxel. Must be between 0 and 1.
    #' @param axis_concentration A numeric value specifying the concentration of
    #'   cylinders along the mean axis. Defaults to `Inf`.
    #' @param radius_sd A numeric value specifying the standard deviation of the
    #'  radius of the cylinder population in meters. Defaults to `0`.
    #' @param radial_model A character string specifying the radial model to
    #'   use. Choices are `"soderman"`, `"callaghan"`, `"stanisz"`, `"neuman"`,
    #'   and `"vangelderen"`. Defaults to `"soderman"`.
    #' @param voxel_size A length-3 numeric vector specifying the size of the
    #'   voxel in meters. Defaults to `c(2e-3, 2e-3, 2e-3)`.
    #'
    #' @return An instance of the [`CylinderBundleCompartment`] class.
    initialize = function(axis,
                          radius,
                          diffusivity,
                          axial_diffusivity,
                          radial_diffusivity,
                          cylinder_density,
                          axis_concentration = Inf,
                          radius_sd = 0,
                          radial_model = c("soderman", "callaghan", "stanisz",
                                           "neuman", "vangelderen"),
                          voxel_size = rep(2e-3, 3)) {
      radial_model <- rlang::arg_match(radial_model)

      private$axis <- axis
      if (cylinder_density < 0 || cylinder_density > 1) {
        cli::cli_abort("The cylinder density must be between 0 and 1.")
      }
      private$cylinder_density <- cylinder_density
      private$axial_diffusivity <- axial_diffusivity
      private$radial_diffusivity <- radial_diffusivity

      voxel_volume <- prod(voxel_size)
      L <- min(voxel_size)
      n_cylinders <- round(voxel_volume * cylinder_density /
        (pi * radius[1]^2 * L), digits = 0)
      cli::cli_alert_info("Number of cylinders: {n_cylinders}")

      withr::with_seed(1234, {
        axis_sample <- rwatson(n_cylinders, axis, axis_concentration)
        radius_sample <- rgamma(n_cylinders, radius, radius_sd)
      })

      private$cylinder_compartments <- purrr::map2(
        .x = axis_sample,
        .y = radius_sample,
        .f = \(axis, radius) {
        CylinderCompartment$new(
          axis = axis,
          radius = radius,
          diffusivity = diffusivity,
          radial_model = radial_model
        )
      })
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
    #' cylinderBundleComp <- CylinderBundleCompartment$new(
    #'   axis = c(0, 0, 1),
    #'   radius = 1e-5,
    #'   diffusivity = 2.0e-9,
    #'   cylinder_density = 0.5,
    #'   axial_diffusivity = 2.0e-9,
    #'   radial_diffusivity = 2.0e-10,
    #'   radial_model = "soderman",
    #'   voxel_size = c(1, 1, 1) * 1e-3
    #' )
    #' cylinderBundleComp$get_signal(
    #'   small_delta = 0.03,
    #'   big_delta = 0.03,
    #'   G = 0.040,
    #'   direction = c(0, 0, 1)
    #' )
    get_signal = function(small_delta, big_delta, G, direction,
                          echo_time = NULL,
                          n_max = 20L,
                          m_max = 50L) {
      cylinder_contribs <- purrr::map_dbl(
        .x = private$cylinder_compartments,
        .f = \(compartment) {
        compartment$get_signal(
          small_delta = small_delta,
          big_delta = big_delta,
          G = G,
          direction = direction,
          echo_time = echo_time,
          n_max = n_max,
          m_max = m_max
        )
      })
      restricted_signal <- mean(cylinder_contribs)
      bvalue <- private$gamma^2 * small_delta^2 * G^2 * (big_delta - small_delta / 3)
      work_value <- (private$axial_diffusivity - private$radial_diffusivity) * sum(direction * private$axis)^2
      hindered_signal <- exp(-bvalue * (private$radial_diffusivity + work_value))
      return(private$cylinder_density * restricted_signal +
               (1 - private$cylinder_density) * hindered_signal)
    }
  ),
  private = list(
    axis = NULL,
    cylinder_compartments = NULL,
    cylinder_density = NULL,
    axial_diffusivity = NULL, # m^2.s^-1
    radial_diffusivity = NULL, # m^2.s^-1
    gamma = 2.675987e8 # rad.s^-1.T^-1
  )
)

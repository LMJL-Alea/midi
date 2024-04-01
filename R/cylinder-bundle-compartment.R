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
    #' @param axis A numeric vector of length 3 and unit norm specifying the
    #'   mean axis of the cylinder population.
    #' @param radius A positive numeric value specifying the mean radius of the
    #'   cylinder population in meters.
    #' @param diffusivity A positive numeric value specifying the diffusivity
    #'   within the cylinders in m\eqn{^2}.s\eqn{^{-1}}.
    #' @param cylinder_density A numeric value specifying the density of the
    #'  cylinders in the voxel. Must be between 0 and 1.
    #' @param axial_diffusivity A numeric value specifying the axial diffusivity
    #'   in the space outside the cylinders in m\eqn{^2}.s\eqn{^{-1}}. If not
    #'   provided, defaults to a tortuosity model reducing the intrinsic
    #'   diffusivity depending on orientation dispersion. Defaults to `NULL`.
    #' @param radial_diffusivity A numeric value specifying the radial
    #'   diffusivity in the space outside the cylinders in
    #'   m\eqn{^2}.s\eqn{^{-1}}. If not provided, defaults to a tortuosity model
    #'   reducing the axial diffusivity depending on radius heterogeneity.
    #'   Defaults to `NULL`.
    #' @param n_cylinders An integer value specifying the number of cylinders in
    #'   the bundle. Defaults to `1L`.
    #' @param axis_concentration A numeric value specifying the concentration of
    #'   cylinders along the mean axis. Defaults to `Inf`.
    #' @param radius_sd A numeric value specifying the standard deviation of the
    #'  radius of the cylinder population in meters. Defaults to `0`.
    #' @param radial_model A character string specifying the radial model to
    #'   use. Choices are `"soderman"`, `"callaghan"`, `"stanisz"`, `"neuman"`,
    #'   and `"vangelderen"`. Defaults to `"soderman"`.
    #' @param seed An integer value specifying the seed for the random number
    #'   generator. Defaults to `1234`.
    #'
    #' @return An instance of the [`CylinderBundleCompartment`] class.
    initialize = function(axis,
                          radius,
                          diffusivity,
                          cylinder_density,
                          axial_diffusivity = NULL,
                          radial_diffusivity = NULL,
                          n_cylinders = 1L,
                          axis_concentration = Inf,
                          radius_sd = 0,
                          radial_model = c("soderman", "callaghan", "stanisz",
                                           "neuman", "vangelderen"),
                          seed = 1234) {
      radial_model <- rlang::arg_match(radial_model)

      # Set axis
      if (length(axis) != 3L || abs(sum(axis^2) - 1) > sqrt(.Machine$double.eps))
        cli::cli_abort("The axis must be a length-3 unit vector.")
      private$axis <- axis

      # Set radius
      if (radius <= 0)
        cli::cli_abort("The radius must be positive.")
      private$radius <- radius

      # Set diffusivity
      if (diffusivity <= 0)
        cli::cli_abort("The diffusivity must be positive.")
      private$diffusivity <- diffusivity

      # Set cylinder density
      if (cylinder_density < 0 || cylinder_density > 1) {
        cli::cli_abort("The cylinder density must be between 0 and 1.")
      }
      private$cylinder_density <- cylinder_density

      # Set axial diffusivity
      if (is.null(axial_diffusivity)) {
        if (is.infinite(axis_concentration)) {
          axial_diffusivity <- diffusivity
        } else {
          axial_diffusivity <- diffusivity * concentration_index(axis_concentration)^2
        }
      } else {
        if (axial_diffusivity > diffusivity)
          cli::cli_abort("The axial diffusivity should not be greater than the intrinsic diffusivity.")
      }
      private$axial_diffusivity <- axial_diffusivity

      # Set radial diffusivity
      if (is.null(radial_diffusivity)) {
        radial_diffusivity <- axial_diffusivity * (1 - 0.8 * cylinder_density)
      } else {
        if (radial_diffusivity > axial_diffusivity)
          cli::cli_abort("The radial diffusivity should not be greater than the axial diffusivity.")
      }
      private$radial_diffusivity <- radial_diffusivity

      if (cylinder_density > 0) {
        if (is.infinite(axis_concentration) && radius_sd == 0) {
          axis_sample <- list(axis)
          radius_sample <- list(radius)
        } else {
          withr::with_seed(seed, {
            axis_sample <- rwatson(n_cylinders, axis, axis_concentration)
            radius_sample <- rgamma(n_cylinders, radius, radius_sd)
          })
        }

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
      }
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
    #'   radial_model = "soderman"
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
      b <- bvalue(small_delta, big_delta, G)
      work_value <- (private$axial_diffusivity - private$radial_diffusivity) * sum(direction * private$axis)^2
      hindered_signal <- exp(-b * (private$radial_diffusivity + work_value))
      if (private$cylinder_density < .Machine$double.eps) {
        return(hindered_signal)
      }
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
      return(private$cylinder_density * restricted_signal +
               (1 - private$cylinder_density) * hindered_signal)
    }
  ),
  private = list(
    axis = NULL,
    radius = NULL,
    diffusivity = NULL,
    cylinder_compartments = NULL,
    cylinder_density = NULL,
    axial_diffusivity = NULL, # m^2.s^-1
    radial_diffusivity = NULL # m^2.s^-1
  )
)

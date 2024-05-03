#' Cylinder bundle compartment class
#'
#' @description A class to model restricted diffusion in a bundle of cylinders.
#'
#' @export
#' @examples
#' withr::with_seed(1234, {
#'   cyls <- rcylinders(
#'     n = 100,
#'     axis_mean = c(0, 0, 1),
#'     radius_mean = 5,
#'     diffusivity_mean = 3,
#'     axis_concentration = 10,
#'     radius_sd = 1,
#'     diffusivity_sd = 0
#'   )
#' })
#'
#' comp <- CylinderBundleCompartment$new(
#'   cylinder_density = 0.5,
#'   cylinder_compartments = cyls
#' )
#'
#' comp$get_signal(
#'   small_delta = 30,
#'   big_delta = 30,
#'   G = 0.040,
#'   direction = c(0, 0, 1)
#' )
#'
#' comp$get_parameters()
CylinderBundleCompartment <- R6::R6Class(
  "CylinderBundleCompartment",
  inherit = BaseCompartment,
  public = list(
    #' @description Instantiates a new cylinder bundle compartment.
    #'
    #' @param cylinder_density A numeric value specifying the density of the
    #'  cylinders in the voxel. Must be between 0 and 1.
    #' @param cylinder_compartments A list of instances of the
    #'   [`CylinderCompartment`] class specifying the compartments within the
    #'   bundle.
    #' @param axial_diffusivity A numeric value specifying the axial diffusivity
    #'   in the space outside the cylinders in m\eqn{^2}.s\eqn{^{-1}}. If not
    #'   provided, defaults to a tortuosity model reducing the intrinsic
    #'   diffusivity depending on orientation dispersion. Defaults to `NULL` in
    #'   which case the extracellular axial diffusivity is set via a tortuosity
    #'   model based on the dispersion in orientation.
    #' @param radial_diffusivity A numeric value specifying the radial
    #'   diffusivity in the space outside the cylinders in
    #'   m\eqn{^2}.s\eqn{^{-1}}. Defaults to `NULL` in which case the
    #'   extracellular radial diffusivity is set via a tortuosity model based on
    #'   the intracellular density.
    #'
    #' @return An instance of the [`CylinderBundleCompartment`] class.
    initialize = function(cylinder_density,
                          cylinder_compartments,
                          axial_diffusivity = NULL,
                          radial_diffusivity = NULL) {
      # Set cylinder density
      if (cylinder_density < 0 || cylinder_density > 1) {
        cli::cli_abort("The cylinder density must be between 0 and 1.")
      }
      private$cylinder_density <- cylinder_density

      # Set cylinder compartments
      private$cylinder_compartments <- cylinder_compartments

      # Retrieve parameters from the cylinder compartments
      params <- cylinder_compartments |>
        purrr::map(\(comp) {
          if (!inherits(comp, "CylinderCompartment")) {
            cli::cli_abort("All compartments must be of class CylinderCompartment.")
          }
          comp$get_parameters()
        }) |>
        purrr::list_transpose()

      # Then, extract axis sample and fit Watson distribution
      wd <- WatsonDistribution$new()
      wd$fit(params$CylinderAxis)
      private$axis <- wd$get_axis()
      private$axis_concentration <- wd$get_concentration()

      gd <- GammaDistribution$new()
      # Then, extract radius sample and fit gamma distribution
      radius_sample <- params$CylinderRadius
      if (var(radius_sample) == 0) {
        private$radius <- radius_sample[1]
        private$radius_sd <- 0
      } else {
        gd$fit(params$CylinderRadius)
        private$radius <- gd$get_shape() * gd$get_scale()
        private$radius_sd <- sqrt(gd$get_shape()) * gd$get_scale()
      }

      # Finally, extract diffusivity sample and fit gamma distribution
      diffusivity_sample <- params$CylinderDiffusivity
      if (var(diffusivity_sample) == 0) {
        private$diffusivity <- diffusivity_sample[1]
        private$diffusivity_sd <- 0
      } else {
        gd$fit(params$CylinderDiffusivity)
        private$diffusivity <- gd$get_shape() * gd$get_scale()
        private$diffusivity_sd <- sqrt(gd$get_shape()) * gd$get_scale()
      }

      # Set axial diffusivity
      if (is.null(axial_diffusivity)) {
        if (is.infinite(private$axis_concentration)) {
          axial_diffusivity <- private$diffusivity
        } else {
          axial_diffusivity <- private$diffusivity * wd$get_rvalue()^2
        }
      } else {
        if (axial_diffusivity > private$diffusivity)
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
    }
  ),
  private = list(
    axis = NULL,
    axis_concentration = NULL,
    radius = NULL,
    radius_sd = NULL,
    diffusivity = NULL,
    diffusivity_sd = NULL,
    cylinder_compartments = NULL,
    cylinder_density = NULL,
    axial_diffusivity = NULL, # um^2.ms^-1
    radial_diffusivity = NULL, # um^2.ms^-1
    compute_signal = function(small_delta, big_delta, G,
                              direction = c(0, 0, 1),
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
      private$cylinder_density * restricted_signal +
        (1 - private$cylinder_density) * hindered_signal
    },
    parameters = function() list(
      CylinderAxisMean = private$axis,
      CylinderAxisConcentration = private$axis_concentration,
      CylinderRadiusMean = private$radius,
      CylinderRadiusStd = private$radius_sd,
      CylinderDiffusivityMean = private$diffusivity,
      CylinderDiffusivityStd = private$diffusivity_sd,
      CylinderDensity = private$cylinder_density,
      CylinderAxialDiffusivity = private$axial_diffusivity,
      CylinderRadialDiffusivity = private$radial_diffusivity
    )
  )
)

#' Sample cylinder compartments
#'
#' This function samples \eqn{n} cylinder compartments with given axis, radius
#' and diffusivity distributions.
#'
#' The axis distribution is given by a mean and a concentration parameter and
#' the Dimroth-Watson distribution is used to sample values. The radius and
#' diffusivity distributions are given by a mean and a standard deviation and
#' the Gamma distribution is used to sample values. If the concentration
#' parameter is set to `Inf` the axis is fixed to the mean value. If the
#' standard deviation of the radius and diffusivity distributions are set to `0`
#' the radius and diffusivity are fixed to the mean values. If all parameters
#' are fixed, only one compartment is sampled.
#'
#' @param n An integer value specifying the number of compartments to sample.
#' @param axis_mean A numeric value specifying the mean of the axis
#'   distribution.
#' @param radius_mean A numeric value specifying the mean of the radius
#'   distribution.
#' @param diffusivity_mean A numeric value specifying the mean of the
#'   diffusivity distribution.
#' @param axis_concentration A numeric value specifying the concentration of the
#'   axis distribution. Defaults to `Inf` in which case the axis is fixed to the
#'   mean value.
#' @param radius_sd A numeric value specifying the standard deviation of the
#'   radius distribution. Defaults to `0` in which case the radius is fixed to
#'   the mean value.
#' @param diffusivity_sd A numeric value specifying the standard deviation of
#'   the diffusivity distribution. Defaults to `0` in which case the diffusivity
#'   is fixed to the mean value.
#'
#' @return A list of \eqn{n} cylinder compartments of class
#'   [`CylinderCompartment`].
#' @export
#'
#' @examples
#' # Sample 10 cylinder compartments with fixed axis, radius and diffusivity
#' # set.seed(42)
#' cyl_distr <- rcylinders(
#'   n = 10L,
#'   axis_mean = c(0, 0, 1),
#'   radius_mean = 5,
#'   diffusivity_mean = 3
#' )
rcylinders <- function(n, axis_mean, radius_mean, diffusivity_mean,
                       axis_concentration = Inf,
                       radius_sd = 0,
                       diffusivity_sd = 0) {
  if (is.infinite(axis_concentration) && radius_sd == 0 && diffusivity_sd == 0) {
    n <- 1L
  }

  # Produce axis sample
  if (is.infinite(axis_concentration)) {
    axis_sample <- purrr::map(1:n, \(.n) axis_mean)
  } else {
    wd <- WatsonDistribution$new(
      mu = axis_mean,
      kappa = axis_concentration
    )
    axis_sample <- wd$random(n)
  }

  # Produce radius sample
  if (radius_sd == 0) {
    radius_sample <- rep(radius_mean, n)
  } else {
    scale <- radius_sd^2 / radius_mean
    shape <- radius_mean / scale
    gd <- GammaDistribution$new(
      shape = shape,
      scale = scale
    )
    radius_sample <- gd$random(n)
  }

  # Produce diffusivity sample
  if (diffusivity_sd == 0) {
    diffusivity_sample <- rep(diffusivity_mean, n)
  } else {
    scale <- diffusivity_sd^2 / diffusivity_mean
    shape <- diffusivity_mean / scale
    gd <- GammaDistribution$new(
      shape = shape,
      scale = scale
    )
    diffusivity_sample <- gd$random(n)
  }

  # Creates cylinder compartments
  purrr::pmap(
    .l = list(axis_sample, radius_sample, diffusivity_sample),
    .f = \(.axis, .radius, .diffusivity) {
    restr_comp <- VanGelderenCompartment$new(
      radius = .radius,
      diffusivity = .diffusivity
    )
    CylinderCompartment$new(
      axis = .axis,
      restricted_compartment = restr_comp
    )
  })
}

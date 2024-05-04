#' @title Base distribution class
#'
#' @description This class defines the base distribution class from which all
#'   other distributions inherit. It provides a common interface for fitting the
#'   distribution to data and generating random samples.
#'
#' @keywords internal
BaseDistribution <- R6::R6Class(
  classname = "BaseDistribution",
  public = list(
    #' @description Fit the distribution to the data.
    #'
    #' @param x A numeric vector of data to fit the distribution to.
    #'
    #' @return Internally sets the parameters of the distribution. Returns
    #'   nothing.
    fit = function(x) {
      private$fit_impl(x)
      invisible(self)
    },

    #' @description Generate random samples from the distribution.
    #'
    #' @param n An integer specifying the number of samples to generate.
    #'
    #' @return A numeric vector of size `n` containing the random samples.
    random = function(n) {
      private$random_impl(n)
    },

    #' @description Compute the mean of the distribution.
    #'
    #' @return A numeric value storing the mean of the distribution.
    get_mean = function() {
      private$mean_impl()
    },

    #' @description Compute the variance of the distribution.
    #'
    #' @return A numeric value storing the variance of the distribution.
    get_variance = function() {
      private$variance_impl()
    }
  )
)

#' Gamma distribution class
#'
#' @description This class defines the gamma distribution. It provides methods
#'   for fitting the distribution to data and generating random samples.
#'
#' @export
#' @examples
#' gd <- GammaDistribution$new(
#'   shape = 1,
#'   scale = 5
#' )
#' gd$get_shape()
#' gd$get_scale()
#' gd$get_mean()
#' gd$get_variance()
#' gd$random(10)
#' gd$fit(gd$random(100))
#' gd$get_shape()
#' gd$get_scale()
GammaDistribution <- R6::R6Class(
  classname = "GammaDistribution",
  inherit = BaseDistribution,
  public = list(
    #' @description Creates a new gamma distribution.
    #'
    #' @param shape A numeric value specifying the shape parameter of the gamma
    #'   distribution. Defaults to `1.0`.
    #' @param scale A numeric value specifying the scale parameter of the gamma
    #'   distribution. Defaults to `1.0`.
    #'
    #' @return An instance of the gamma distribution as an object of class
    #'   [`GammaDistribution`].
    initialize = function(shape = 1.0, scale = 1.0) {
      private$shape <- shape
      private$scale <- scale
    },

    #' @description Retrieves the shape parameter of the gamma distribution.
    #'
    #' @return A numeric value storing the shape parameter of the gamma
    #'   distribution.
    get_shape = function() {
      private$shape
    },

    #' @description Retrieves the scale parameter of the gamma distribution.
    #'
    #' @return A numeric value storing the scale parameter of the gamma
    #'   distribution.
    get_scale = function() {
      private$scale
    }
  ),
  private = list(
    shape = NULL,
    scale = NULL,
    fit_impl = function(x) {
      n <- length(x)
      lnx <- log(x)
      xlnx <- x * lnx
      x_bar <- mean(x)
      lnx_bar <- mean(lnx)
      xlnx_bar <- mean(xlnx)
      theta <- xlnx_bar - x_bar * lnx_bar
      kappa <- x_bar / theta
      private$shape <- kappa - (3 * kappa - 2 * kappa / (1 + kappa) / 3 -
                                  4 * kappa / (1 + kappa)^2 / 5) / n
      private$scale <- n * theta / (n - 1)
    },
    random_impl = function(n) {
      if (private$scale < .Machine$double.eps) return(rep(mean(), n))
      stats::rgamma(n, shape = private$shape, scale = private$scale)
    },
    mean_impl = function() {
      private$shape * private$scale
    },
    variance_impl = function() {
      private$shape * private$scale^2
    }
  )
)

#' Watson distribution class
#'
#' @description This class defines the Watson distribution. It provides methods
#'   for fitting the distribution to data and generating random samples.
#'
#' @export
#' @examples
#' wd <- WatsonDistribution$new(
#'   mu = c(0, 0, 1),
#'   kappa = 10
#' )
#' wd$get_axis()
#' wd$get_concentration()
#' wd$random(10)
#' wd$fit(wd$random(100))
#' wd$get_axis()
#' wd$get_concentration()
WatsonDistribution <- R6::R6Class(
  classname = "WatsonDistribution",
  inherit = BaseDistribution,
  public = list(
    #' @description Creates a new Watson distribution.
    #'
    #' @param mu A numeric vector of length 3 specifying the mean direction of
    #'   the Watson distribution. Defaults to `(0, 0, 1)`.
    #' @param kappa A numeric value specifying the concentration parameter of
    #'   the Watson distribution. Defaults to `10.0`.
    #'
    #' @return An instance of the Watson distribution as an object of class
    #'   [`WatsonDistribution`].
    initialize = function(mu = c(0, 0, 1), kappa = 10.0) {
      private$mu <- mu
      private$kappa <- kappa
    },

    #' @description Retrieves the mean axis of the Watson distribution.
    #'
    #' @return A numeric vector of length 3 storing the mean axis of the Watson
    #'   distribution.
    get_axis = function() {
      private$mu
    },

    #' @description Retrieves the concentration parameter of the Watson
    #'   distribution.
    #'
    #' @return A numeric value storing the concentration parameter of the Watson
    #'   distribution.
    get_concentration = function() {
      private$kappa
    }
  ),
  private = list(
    mu = NULL,
    kappa = NULL,
    r = 1.0,
    fit_impl = function(x) {
      n <- length(x)
      x <- do.call(rbind, x)
      dcm <- t(x) %*% x
      eig <- eigen(dcm)
      private$mu <- eig$vectors[, which.max(eig$values)]
      private$r <- max(eig$values) / n
      if (abs(private$r - 1) < sqrt(.Machine$double.eps)) {
        private$kappa <- Inf
      } else {
        private$kappa <- stats::uniroot(
          f = \(.x) concentration_index(.x) - private$r,
          interval = c(0, 100)
        )$root
      }
    },
    random_impl = function(n) {
      if (private$kappa == Inf) return(purrr::map(1:n, \(.n) private$mu))

      # Computes rotation matrix that brings z-axis to mu
      R <- rotation_matrix_from_z_to_mu(private$mu)

      purrr::map(1:n, \(.n) {
        if (private$kappa > sqrt(.Machine$double.eps)) {
          U <- stats::runif(1)
          S <- 1 + log(U + (1 - U) * exp(-private$kappa)) / private$kappa
          V <- stats::runif(1)
          if (V > 1e-6) {
            while (log(V) > private$kappa * S * (S - 1)) {
              U <- stats::runif(1)
              S <- 1 + log(U + (1 - U) * exp(-private$kappa)) / private$kappa
              V <- stats::runif(1)
              if (V < 1e-6) {
                break
              }
            }
          }
        } else if (private$kappa < -sqrt(.Machine$double.eps)) {
          C1 <- sqrt(abs(private$kappa))
          C2 <- atan(C1)
          U <- stats::runif(1)
          V <- stats::runif(1)
          S <- (1 / C1) * tan(C2 * U)
          T <- private$kappa * S * S
          while (V > (1 - T) * exp(T)) {
            U <- stats::runif(1)
            V <- stats::runif(1)
            S <- (1 / C1) * tan(C2 * U)
            T <- private$kappa * S * S
          }
        } else {
          S <- cos(pi * stats::runif(1))
        }

        phi <- 2 * pi * stats::runif(1)

        out <- c(sqrt(1 - S * S) * cos(phi), sqrt(1 - S * S) * sin(phi), S)
        out <- R %*% out
        out <- out / sqrt(sum(out^2))
        as.numeric(out)
      })
    },
    mean_impl = function() {
      private$mu
    },
    variance_impl = function() {
      1 - private$r
    }
  )
)

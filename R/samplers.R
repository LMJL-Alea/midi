# Computes rotation matrix that brings z-axis to mu
rotation_matrix_from_z_to_mu <- function(mu) {
  mu <- mu / sqrt(sum(mu^2))
  z <- c(0, 0, 1)
  if (all(mu == z)) {
    return(diag(3))
  }
  if (all(mu == -z)) {
    return(matrix(c(1, 0, 0, 0, -1, 0, 0, 0, -1), nrow = 3L, ncol = 3L))
  }
  # computes cross product between z and mu
  v <- c(
    z[2] * mu[3] - z[3] * mu[2],
    z[3] * mu[1] - z[1] * mu[3],
    z[1] * mu[2] - z[2] * mu[1]
  )
  # computes norm of v
  s <- sqrt(sum(v^2))
  # computes inner product between z and mu
  c <- sum(z * mu)

  V <- matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), nrow = 3, ncol = 3)
  diag(3) + V + V %*% V * (1 - c) / (s^2)
}

rwatson <- function(n, mu, kappa) {
  if (kappa == Inf) return(purrr::map(1:n, \(.n) mu))

  # Computes rotation matrix that brings z-axis to mu
  R <- rotation_matrix_from_z_to_mu(mu)

  purrr::map(1:n, \(.n) {
    if (kappa > sqrt(.Machine$double.eps)) {
      U <- stats::runif(1)
      S <- 1 + log(U + (1 - U) * exp(-kappa)) / kappa
      V <- stats::runif(1)
      if (V > 1e-6) {
        while (log(V) > kappa * S * (S - 1)) {
          U <- stats::runif(1)
          S <- 1 + log(U + (1 - U) * exp(-kappa)) / kappa
          V <- stats::runif(1)
          if (V < 1e-6) {
            break
          }
        }
      }
    } else if (kappa < -sqrt(.Machine$double.eps)) {
      C1 <- sqrt(abs(kappa))
      C2 <- atan(C1)
      U <- stats::runif(1)
      V <- stats::runif(1)
      S <- (1 / C1) * tan(C2 * U)
      T <- kappa * S * S
      while (V > (1 - T) * exp(T)) {
        U <- stats::runif(1)
        V <- stats::runif(1)
        S <- (1 / C1) * tan(C2 * U)
        T <- kappa * S * S
      }
    } else {
      S <- cos(pi * stats::runif(1))
    }

    phi <- 2 * pi * stats::runif(1)

    out <- c(sqrt(1 - S * S) * cos(phi), sqrt(1 - S * S) * sin(phi), S)
    out <- R %*% out
    out <- out / sqrt(sum(out^2))
    out
  })
}

rgamma <- function(n, mu, sd) {
  if (sd < .Machine$double.eps) return(rep(mu, n))
  shape <- mu^2 / sd^2
  scale <- sd^2 / mu
  stats::rgamma(n, shape = shape, scale = scale)
}

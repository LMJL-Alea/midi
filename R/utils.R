bvalue <- function(small_delta, big_delta, G) {
  gyromagnetic_ratio()^2 * small_delta^2 * G^2 * (big_delta - small_delta / 3)
}

besselJ_derivative <- function(x, n) {
  ifelse(
    n == 0,
    -besselJ(x, 1),
    (besselJ(x, n - 1) - besselJ(x, n + 1)) / 2
  )
}

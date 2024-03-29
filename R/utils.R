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

# this is 2 exp(-x) M(x, 0.5, 1.5)
kummer <- Vectorize(function(x) {
  integrate(\(t) exp(-x * (1 - t)) * t^-0.5, 0, 1)$value
})
# this is 2 exp(-x) M'(x, 0.5, 1.5)
kummer_derivative <- Vectorize(function(x) {
  integrate(\(t) exp(-x * (1 - t)) * t^0.5, 0, 1)$value
})

concentration_index <- function(k) {
  kummer_derivative(k) / kummer(k)
}

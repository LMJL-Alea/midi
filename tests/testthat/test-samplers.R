test_that("`rotation_matrix_from_z_to_mu()` works", {
  z <- c(0, 0, 1)
  R <- rotation_matrix_from_z_to_mu(c(0, 0, 1))
  expect_equal(as.numeric(R %*% z), z)
  R <- rotation_matrix_from_z_to_mu(c(1, 0, 0))
  expect_equal(as.numeric(R %*% z), c(1, 0, 0))
  R <- rotation_matrix_from_z_to_mu(c(0, 0, -1))
  expect_equal(as.numeric(R %*% z), c(0, 0, -1))
})

test_that("`rwatson()` works", {
  mu <- c(0, 0, 1)
  n <- 10L
  kappa <- Inf
  samples <- rwatson(n, mu, kappa)
  expect_equal(length(samples), n)
  kappa <- 10
  samples <- rwatson(n, mu, kappa)
  expect_equal(length(samples), n)
  kappa <- -10
  samples <- rwatson(n, mu, kappa)
  expect_equal(length(samples), n)
  kappa <- 0
  samples <- rwatson(n, mu, kappa)
  expect_equal(length(samples), n)
})

test_that("`rgamma() works", {
  mu <- 1
  n <- 10L
  sd <- 1
  samples <- rgamma(n, mu, sd)
  expect_equal(length(samples), n)
  sd <- 0
  samples <- rgamma(n, mu, sd)
  expect_equal(length(samples), n)
})

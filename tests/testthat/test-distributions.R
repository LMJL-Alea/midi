test_that("`rotation_matrix_from_z_to_mu()` works", {
  z <- c(0, 0, 1)
  R <- rotation_matrix_from_z_to_mu(c(0, 0, 1))
  expect_equal(as.numeric(R %*% z), z)
  R <- rotation_matrix_from_z_to_mu(c(1, 0, 0))
  expect_equal(as.numeric(R %*% z), c(1, 0, 0))
  R <- rotation_matrix_from_z_to_mu(c(0, 0, -1))
  expect_equal(as.numeric(R %*% z), c(0, 0, -1))
})

test_that("the Watson distribution` works", {
  n <- 10L
  mu <- c(0, 0, 1)
  kappa <- Inf
  wd <- WatsonDistribution$new(mu, kappa)
  samples <- wd$random(n)
  expect_equal(length(samples), n)
  kappa <- 10
  wd <- WatsonDistribution$new(mu, kappa)
  samples <- wd$random(n)
  expect_equal(length(samples), n)
  kappa <- -10
  wd <- WatsonDistribution$new(mu, kappa)
  samples <- wd$random(n)
  expect_equal(length(samples), n)
  kappa <- 0
  wd <- WatsonDistribution$new(mu, kappa)
  samples <- wd$random(n)
  expect_equal(length(samples), n)
})

test_that("the Gamma distribution works", {
  n <- 10L
  shape <- 1
  scale <- 1
  gd <- GammaDistribution$new(shape, scale)
  samples <- gd$random(n)
  expect_equal(length(samples), n)
})

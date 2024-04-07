test_that("`simulate_bundle()` works", {
  density <- 0.5
  voxel_size <- 5 # micrometers
  withr::with_seed(1234, {
    out <- simulate_bundle(density, voxel_size)
  })
  expect_equal(colnames(out$sections), c("x", "y", "r"))
})

test_that("`autplot.bundle()` works", {
  density <- 0.5
  voxel_size <- 5 # micrometers
  withr::with_seed(1234, {
    out <- simulate_bundle(density, voxel_size)
  })
  p <- ggplot2::autoplot(out)
  expect_true(inherits(p, "gg"))
})

test_that("`plot3d()` works", {
  density <- 0.5
  voxel_size <- 5 # micrometers
  withr::with_seed(1234, {
    out <- simulate_bundle(density, voxel_size)
  })
  p <- plot3d(out)
  expect_true(inherits(p, "plotly"))
})

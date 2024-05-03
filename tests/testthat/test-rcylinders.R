test_that("`rcylinders()` works", {
  withr::with_seed(1234, {
    cyl_distr <- rcylinders(
      n = 10L,
      axis_mean = c(0, 0, 1),
      radius_mean = 5,
      diffusivity_mean = 3
    )
  })
  expect_equal(length(cyl_distr), 1L)

  withr::with_seed(1234, {
    cyl_distr <- rcylinders(
      n = 10L,
      axis_mean = c(0, 0, 1),
      axis_concentration = 10,
      radius_mean = 5,
      diffusivity_mean = 3
    )
  })
  expect_equal(length(cyl_distr), 10L)
})

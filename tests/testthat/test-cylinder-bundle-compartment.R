test_that("`CylinderBundleCompartment` class works", {
  withr::with_seed(1234, {
    cyls <- rcylinders(
      n = 100,
      axis_mean = c(0, 0, 1),
      radius_mean = 5,
      diffusivity_mean = 3,
      axis_concentration = 10,
      radius_sd = 1,
      diffusivity_sd = 0
    )
  })

  comp <- CylinderBundleCompartment$new(
    cylinder_density = 0.5,
    cylinder_compartments = cyls
  )

  out <- comp$get_signal(
    small_delta = 30,
    big_delta = 30,
    G = 0.040,
    direction = c(0, 0, 1)
  )

  expect_equal(round(out, 5L), 0.00551)

  params <- comp$get_parameters()

  expect_equal(round(params$CylinderAxisMean, 5L), c(-0.02893, -0.00652, -0.99956))
  expect_equal(round(params$CylinderAxisConcentration, 5L), 10.68640)
  expect_equal(round(params$CylinderRadiusMean, 5L), 5.08216)
  expect_equal(round(params$CylinderRadiusStd, 5L), 0.99304)
  expect_equal(round(params$CylinderDiffusivityMean, 5L), 3.00000)
  expect_equal(round(params$CylinderDiffusivityStd, 5L), 0.00000)
  expect_equal(params$CylinderDensity, 0.5)
  expect_equal(round(params$CylinderAxialDiffusivity, 5L), 2.43153)
  expect_equal(round(params$CylinderRadialDiffusivity, 5L), 1.45892)
})

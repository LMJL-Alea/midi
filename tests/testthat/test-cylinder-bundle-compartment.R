test_that("`CylinderBundleCompartment` class works", {
  cbComp <- CylinderBundleCompartment$new(
    axis = c(0, 0, 1),
    radius = 1e-5,
    diffusivity = 2.0e-9,
    cylinder_density = 0.5,
    radial_model = "soderman"
  )
  out <- cbComp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040,
    direction = c(0, 0, 1)
  )
  expect_true(inherits(cbComp, "CylinderBundleCompartment"))
})

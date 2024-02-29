test_that("`CylinderCompartment` class works", {
  cylinderComp <- CylinderCompartment$new(
    axis = c(0, 0, 1),
    radius = 1e-6,
    diffusivity = 2.0e-9,
    radial_model = "soderman"
  )
  out <- cylinderComp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040,
    direction = c(0, 0, 1)
  )
  expect_true(inherits(cylinderComp, "CylinderCompartment"))
})

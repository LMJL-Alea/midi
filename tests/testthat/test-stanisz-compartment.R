test_that("`StaniszCompartment` class works", {
  comp <- StaniszCompartment$new(
    radius = 1e-6,
    diffusivity = 2.0e-9
  )
  out <- comp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040
  )
  expect_true(inherits(comp, "CylinderRadialCompartment"))
})

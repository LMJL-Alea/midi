test_that("`NeumanCompartment` class works", {
  comp <- NeumanCompartment$new(
    radius = 1e-6,
    diffusivity = 2.0e-9
  )
  out <- comp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040,
    echo_time = 0.080
  )
  expect_true(inherits(comp, "CylinderRadialCompartment"))
})

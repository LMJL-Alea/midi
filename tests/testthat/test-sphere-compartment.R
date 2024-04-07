test_that("`SphereCompartment` class works", {
  sphereComp <- SphereCompartment$new(
    radius = 1e-6,
    diffusivity = 2.0e-9,
    model = "soderman"
  )
  out <- sphereComp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040
  )
  expect_true(inherits(sphereComp, "SphereCompartment"))
})

test_that("`SphereCompartment` class works", {
  sphereComp <- SphereCompartment$new()
  out <- sphereComp$get_signal(
    small_delta = 30,
    big_delta = 30,
    G = 0.040e-3
  )
  expect_true(inherits(sphereComp, "SphereCompartment"))
  expect_equal(round(out, 5), 0.12045)
})

test_that("`CallaghanCompartment` class works", {
  comp <- CallaghanCompartment$new(
    radius = 5,
    diffusivity = 3
  )
  expect_equal(comp$get_parameters(), list(
    Radius = 5,
    Diffusivity = 3
  ))
  out <- comp$get_signal(
    small_delta = 30,
    big_delta = 30,
    G = 0.040
  )
  expect_true(inherits(comp, "CircularlyShapedCompartment"))
  expect_equal(round(out, 5), 0.50512)
})

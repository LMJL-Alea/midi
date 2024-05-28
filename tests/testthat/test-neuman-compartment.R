test_that("`NeumanCompartment` class works", {
  comp <- NeumanCompartment$new(
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
    G = 0.040,
    echo_time = 80
  )
  expect_true(inherits(comp, "CircularlyShapedCompartment"))
  expect_equal(round(out, 5), 0.93140)
})

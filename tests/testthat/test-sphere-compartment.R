test_that("`SphereCompartment` class works", {
  restr_comp <- VanGelderenCompartment$new(
    radius = 5,
    diffusivity = 3
  )
  comp <- SphereCompartment$new(restr_comp)
  expect_equal(comp$get_parameters(), list(
    SphereRadius = 5,
    SphereDiffusivity = 3
  ))
  out <- comp$get_signal(
    small_delta = 30,
    big_delta = 30,
    G = 0.040
  )
  expect_true(inherits(comp, "CircularlyShapedCompartment"))
  expect_equal(round(out, 5), 0.91246)
})

test_that("`CylinderCompartment` class works", {
  restr_comp <- VanGelderenCompartment$new(
    radius = 5,
    diffusivity = 3
  )
  comp <- CylinderCompartment$new(
    axis = c(0, 0, 1),
    restricted_compartment = restr_comp
  )
  expect_equal(comp$get_parameters(), list(
    CylinderAxis = c(0, 0, 1),
    CylinderRadius = 5,
    CylinderDiffusivity = 3
  ))
  out <- comp$get_signal(
    small_delta = 30,
    big_delta = 30,
    G = 0.040,
    direction = c(1, 0, 0)
  )
  expect_true(inherits(comp, "CylinderCompartment"))
  expect_equal(round(out, 5), 0.91246)
})

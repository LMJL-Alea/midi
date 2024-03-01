test_that("`CylinderCompartment` class works", {
  cbComp <- CylinderBundleCompartment$new(
    axis = c(0, 0, 1),
    radius = 1e-5,
    diffusivity = 2.0e-9,
    cylinder_density = 0.5,
    axial_diffusivity = 2.0e-9,
    radial_diffusivity = 2.0e-10,
    radial_model = "soderman",
    voxel_size = c(2, 2, 2) * 1e-3
  )
  out <- cbComp$get_signal(
    small_delta = 0.03,
    big_delta = 0.03,
    G = 0.040,
    direction = c(0, 0, 1)
  )
  expect_true(inherits(cbComp, "CylinderBundleCompartment"))
})

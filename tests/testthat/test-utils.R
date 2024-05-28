test_that("bvalue() works", {
  actual <- round(bvalue(small_delta = 13, big_delta = 22, G = 0.040), 3)
  expected <- 0.342
  expect_equal(actual, expected)
})

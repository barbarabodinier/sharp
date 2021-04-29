test_that("output of SimulateGraphical() is of correct dimension", {
  expect_equal(nrow(SimulateGraphical(n = 100)$data), 100)
})

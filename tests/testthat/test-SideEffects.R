test_that("potential side effects", {
  # Make sure seed is not disturbed
  set.seed(1)
  simul <- SimulateGraphical(n = 50, pk = 20, nu = 0.1)
  set.seed(1)
  a <- runif(1)
  set.seed(1)
  stab <- GraphicalModel(xdata = simul$data, K = 5, Lambda_cardinal = 5, verbose = FALSE)
  expect_equal(runif(1), a)
})

test_that("Clustering() models are working", {
  # Data simulation
  set.seed(1)
  simul <- SimulateGraphical(pk = 50, n = 10)

  stab <- Clustering(xdata = simul$data)
  membership <- Clusters(stab)
  expect_equal(length(membership), nrow(simul$data))
})

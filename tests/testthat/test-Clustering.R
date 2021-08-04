test_that("Clustering() models are working", {
  # Data simulation
  set.seed(1)
  simul <- SimulateGraphical(pk = 5, n = 10)

  stab <- Clustering(xdata = simul$data, Lambda = 1.1, K = 5)
  membership <- Clusters(stab)
  expect_equal(length(membership), nrow(simul$data))
})

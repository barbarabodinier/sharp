test_that("BiSelection() models are working", {
  # Data simulation
  K <- 5
  pk <- 15
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
  colnames(ydata) <- paste0("outcome", 1:4)
  x <- simul$X
  y <- ydata

  # sPLS: sparsity on both X and Y
  stab <- BiSelection(
    xdata = x, ydata = y,
    family = "gaussian", K = K, ncomp = 2,
    LambdaX = 1:2,
    LambdaY = 1:2,
    implementation = "SparsePLS"
  )

  # sgPLS: sparsity on both X and Y
  stab <- BiSelection(
    xdata = x, ydata = y,
    group_x = c(10, 5), group_y = c(1, 3),
    family = "gaussian", K = K, ncomp = 2,
    LambdaX = 1:2, AlphaX = c(0.1, 0.3),
    LambdaY = 1:2, AlphaY = c(0.1, 0.3),
    implementation = "SparseGroupPLS"
  )
})

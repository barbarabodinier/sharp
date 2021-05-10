test_that("sgPLS models are working", {
  K <- 5
  pk <- 15
  # Sparse PLS (1 outcome)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
  colnames(ydata) <- paste0("outcome", 1:4)
  x <- simul$X
  y <- ydata

  # With fixed sparsity within groups
  alpha <- 0.5
  stab <- VariableSelection(
    xdata = x, ydata = y, K = K,
    group_x = c(10, 5), alpha.x = alpha,
    Lambda = 1:2,
    implementation = SparseGroupPLS, family = "gaussian",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)

  # Calibrating number of components and groups/sparsity in outcomes
  stab <- BiSelection(
    xdata = x, ydata = y, K = K,
    group_x = c(10, 5),
    group_y = c(1, 3),
    AlphaX = seq(0.1, 0.9, by = 0.1),
    LambdaX = 1:2,
    AlphaY = seq(0.1, 0.9, by = 0.1),
    LambdaY = 1:2,
    ncomp = 3,
    implementation = SparseGroupPLS, family = "gaussian",
    verbose = TRUE
  )
})

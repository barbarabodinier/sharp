test_that("sPLS models are working", {
  K <- 5
  pk <- 7
  # Sparse PLS (1 outcome)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  stab <- VariableSelection(
    xdata = simul$X, ydata = simul$Y, K = K,
    Lambda = 1:(ncol(simul$X) - 1),
    implementation = "SparsePLS", family = "gaussian",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)

  # Sparse PLS (sparse on X, 1 component, multiple outcomes)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
  colnames(ydata) <- paste0("outcome", 1:4)
  x <- simul$X
  y <- ydata
  stab <- VariableSelection(
    xdata = x, ydata = y, K = K,
    Lambda = 1:(ncol(simul$X) - 1), ncomp = 1,
    implementation = "SparsePLS", family = "gaussian",
    verbose = FALSE
  )
  expect_equal(length(SelectedVariables(stab)), pk)
})


test_that("sPLSDA models are working", {
  K <- 5
  pk <- 7

  # Sparse PLS (1 component)
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = pk, family = "binomial")
  Y <- simul$Y
  Y[Y == 0] <- sample(c(0, 2), size = sum(Y == 0), replace = TRUE)
  stab <- VariableSelection(
    xdata = simul$X, ydata = Y, K = K,
    Lambda = 1:(ncol(simul$X) - 1),
    implementation = "SparsePLS", family = "binomial",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)
})

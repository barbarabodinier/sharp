test_that("sPLS models are working", {
  K <- 5
  pk <- 7
  # Sparse PLS (1 outcome)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  stab <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata, K = K,
    Lambda = 1:(ncol(simul$xdata) - 1),
    implementation = SparsePLS, family = "gaussian",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)

  # Sparse PLS (sparse on X, 1 component, multiple outcomes)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = c(3, 2, 2), family = "gaussian")
  x <- simul$xdata
  y <- simul$ydata
  stab <- VariableSelection(
    xdata = x, ydata = y, K = K,
    Lambda = 1:(ncol(simul$xdata) - 1), ncomp = 1,
    implementation = SparsePLS, family = "gaussian",
    verbose = FALSE
  )
  expect_equal(length(SelectedVariables(stab)), pk)
})


test_that("sPLSDA models are working", {
  K <- 5
  pk <- 7

  # Sparse PLS (1 component)
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = c(3, 2, 2), family = "binomial")
  stab <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata, K = K,
    Lambda = 1:(ncol(simul$xdata) - 1),
    implementation = SparsePLS, family = "binomial",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)
})

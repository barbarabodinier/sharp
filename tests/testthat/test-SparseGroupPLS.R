test_that("sgPLS models are working", {
  K <- 10
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
    implementation = "SparseGroupPLS", family = "gaussian",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)

  # Calibrating sparsity within groups
  alpha_list <- seq(0.1, 0.9, by = 0.1)
  selected <- stability_score <- NULL
  for (alpha in alpha_list) {
    stab <- VariableSelection(
      xdata = x, ydata = y, K = K,
      group_x = c(10, 5), alpha.x = alpha,
      Lambda = 1:2,
      implementation = "SparseGroupPLS", family = "gaussian",
      verbose = FALSE
    )
    selected <- rbind(selected, SelectedVariables(stab))
    stability_score <- c(stability_score, max(stab$S, na.rm = TRUE))
  }

  # Calibrating number of components
  stab <- BiSelection(
    xdata = x, ydata = y, K = K,
    group_x = c(10, 5),
    AlphaX = seq(0.1, 0.9, by = 0.1),
    LambdaX = 1:2, ncomp = 3,
    implementation = "SparseGroupPLS", family = "gaussian",
    verbose = TRUE
  )

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
    implementation = "SparseGroupPLS", family = "gaussian",
    verbose = TRUE
  )
})


test_that("sPLSDA models are working", {
  K <- 10
  pk <- 7
  # Sparse PLS-DA (1 outcome)
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "binomial")
  stab <- VariableSelection(
    xdata = simul$X, ydata = simul$Y, K = K,
    Lambda = 1:(ncol(simul$X) - 1),
    implementation = "SparsePLS", family = "gaussian",
    verbose = FALSE
  )
  CalibrationPlot(stab, xlab = "")
  expect_equal(length(SelectedVariables(stab)), pk)

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

  # Sparse PLS (sparse on X, multiple components)
  nvar_x <- selected <- stability_score <- NULL
  for (comp in 1:length(unique(Y))) {
    stab <- VariableSelection(
      xdata = simul$X, ydata = Y, K = K,
      Lambda = 1:(ncol(simul$X) - 1), ncomp = comp,
      implementation = "SparsePLS", family = "binomial",
      verbose = FALSE
    )
    nvar_x <- c(nvar_x, Argmax(stab)[1])
    selected <- rbind(selected, SelectedVariables(stab))
    stability_score <- c(stability_score, max(stab$S, na.rm = TRUE))
  }
  expect_equal(dim(selected), c(3, pk))
  # selected[1:which.max(stability_score), ] # selected in calibrated model
})

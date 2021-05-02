test_that("sPLS models are working", {
  K <- 10
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

  # Sparse PLS (sparse on X, 1 component)
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

  # Sparse PLS (sparse on X, multiple components)
  nvar_x <- selected <- stability_score <- NULL
  for (comp in 1:ncol(y)) {
    stab <- VariableSelection(
      xdata = x, ydata = y, K = K,
      Lambda = 1:(ncol(simul$X) - 1), ncomp = comp,
      keepX_previous = nvar_x,
      implementation = "SparsePLS", family = "gaussian",
      verbose = FALSE
    )
    nvar_x <- c(nvar_x, Argmax(stab)[1])
    selected <- rbind(selected, SelectedVariables(stab))
    stability_score <- c(stability_score, max(stab$S, na.rm = TRUE))
  }
  expect_equal(dim(selected), c(4, pk))
  # selected[1:which.max(stability_score), ] # selected in calibrated model

  # Sparse PLS (sparse on X and Y, multiple components)
  nvar_x <- nvar_y <- selected_x <- selected_y <- stability_score <- NULL
  for (comp in 1:ncol(y)) {
    tmp_stability_score <- tmp_tokeepx <- NULL
    tmp_selected_x <- tmp_selected_y <- NULL
    for (ny in 1:ncol(y)) {
      stab <- VariableSelection(
        xdata = x, ydata = y, K = K,
        Lambda = 1:(ncol(simul$X) - 1), ncomp = comp,
        keepX_previous = nvar_x, keepY = c(nvar_y, ny),
        implementation = "SparsePLS", family = "gaussian",
        verbose = FALSE
      )
      selprop <- apply(Coefficients(stab, side = "Y")[ArgmaxId(stab)[1], , ], 1,
        FUN = function(z) {
          sum(z != 0) / length(z)
        }
      )
      if (any(selprop != 1)) {
        hat_pi <- stab$params$pi_list[which.max(StabilityScore(selprop,
          pi = stab$params$pi_list, K = K
        ))]
        tmp_selected_y <- rbind(
          tmp_selected_y,
          ifelse(selprop >= hat_pi, yes = 1, no = 0)
        )
      } else {
        tmp_selected_y <- rbind(tmp_selected_y, rep(1, ncol(y)))
      }
      tmp_selected_x <- rbind(tmp_selected_x, SelectedVariables(stab))
      tmp_tokeepx <- c(tmp_tokeepx, Argmax(stab)[1])
      tmp_stability_score <- c(tmp_stability_score, max(stab$S, na.rm = TRUE))
    }
    tokeepx <- tmp_tokeepx[which.max(tmp_stability_score)]
    tokeepy <- which.max(tmp_stability_score)
    nvar_x <- c(nvar_x, tokeepx)
    nvar_y <- c(nvar_y, tokeepy)
    selected_x <- rbind(selected_x, tmp_selected_x[which.max(tmp_stability_score), ])
    selected_y <- rbind(selected_y, tmp_selected_y[which.max(tmp_stability_score), ])
    stability_score <- c(stability_score, max(stab$S, na.rm = TRUE))
  }
  # print(selected_x[1:which.max(stability_score), ]) # selected in X in calibrated model
  # print(selected_y[1:which.max(stability_score), ]) # selected in Y in calibrated model
  expect_equal(dim(selected_x), c(4, pk))
  expect_equal(dim(selected_y), c(4, 4))
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
      Lambda = 1:(ncol(simul$X) - 1), ncomp=comp,
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

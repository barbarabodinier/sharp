test_that("lambda grid size for VariableSelection()", {
  skip_on_cran()
  n <- 78
  pk <- 12
  nlambda <- 1
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)
  simul <- SimulateRegression(n = n, pk = pk, family = "gaussian")

  expect_warning(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE
  ))

  stab <- VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda = 0.5,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE
  )
})


test_that("lambda grid size for GraphicalModel()", {
  skip_on_cran()
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- 12
  K <- 5
  nlambda <- 3
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Data simulation
  simul <- SimulateGraphical(n = n, pk = pk)

  expect_warning(GraphicalModel(
    data = simul$data,
    Lambda_cardinal = 1, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE
  ))

  stab <- GraphicalModel(
    data = simul$data,
    Lambda_cardinal = 1, K = K,
    Lambda = 0.5,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE
  )
})

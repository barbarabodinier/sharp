test_that("single-block GraphicalModel()", {
  skip_on_cran()
  for (n_cores in c(1, 2)) {
    PFER_thr <- FDP_thr <- Inf
    n <- 78
    pk <- 12
    K <- 10
    nlambda <- 3
    tau <- 0.55
    n_cat <- 3
    pi_list <- seq(0.6, 0.7, length.out = 15)

    # Data simulation
    set.seed(1)
    simul <- SimulateGraphical(n = n, pk = pk)

    stab <- GraphicalModel(
      xdata = simul$data,
      Lambda_cardinal = nlambda, K = K,
      pi_list = pi_list,
      tau = tau, n_cat = n_cat,
      PFER_thr = PFER_thr,
      FDP_thr = FDP_thr,
      verbose = FALSE,
      start = "cold",
      n_cores = n_cores
    )

    ### Checking dimensions of the outputs
    # Group of outputs 1
    expect_equal(dim(stab$S), c(nlambda, 1))
    expect_equal(dim(stab$Lambda), c(nlambda, 1))
    expect_equal(dim(stab$Q), c(nlambda, 1))
    expect_equal(dim(stab$Q_s), c(nlambda, 1))
    expect_equal(dim(stab$P), c(nlambda, 1))
    expect_equal(dim(stab$PFER), c(nlambda, 1))
    expect_equal(dim(stab$FDP), c(nlambda, 1))

    # Group of outputs 2
    expect_equal(dim(stab$S_2d), c(nlambda, length(pi_list)))
    expect_equal(dim(stab$PFER_2d), c(nlambda, length(pi_list)))
    expect_equal(dim(stab$FDP_2d), c(nlambda, length(pi_list)))

    # Group of outputs 3
    expect_equal(dim(stab$selprop), c(pk, pk, nlambda))
    expect_equal(dim(stab$sign), c(pk, pk))

    # Group of outputs 4
    default_params <- c("graphical_model", "PenalisedGraphical", "cold", "subsampling", FALSE, "MB")
    names(default_params) <- c("type", "implementation", "start", "resampling", "cpss", "PFER_method")
    expect_equal(unlist(stab$methods), default_params)

    # Group of outputs 5
    expect_equal(stab$params$K, K)
    expect_equal(stab$params$pi_list, pi_list)
    expect_equal(stab$params$tau, tau)
    expect_equal(stab$params$n_cat, n_cat)
    expect_equal(stab$params$pk, pk)
    expect_equal(stab$params$PFER_thr, PFER_thr)
    expect_equal(stab$params$FDP_thr, FDP_thr)
    expect_equal(stab$params$n, nrow(simul$data))
  }
})


test_that("multi-block GraphicalModel()", {
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- c(5, 5)
  nlambda <- 4
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Data simulation
  set.seed(1)
  simul <- SimulateGraphical(n = n, pk = pk)

  stab <- GraphicalModel(
    xdata = simul$data, pk = pk,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    start = "cold",
    verbose = FALSE
  )

  ### Checking dimensions of the outputs
  # Group of outputs 1
  expect_equal(dim(stab$S), c(nlambda * 3, 3))
  expect_equal(dim(stab$Lambda), c(nlambda * 3, 3))
  expect_equal(dim(stab$Q), c(nlambda * 3, 3))
  expect_equal(dim(stab$Q_s), c(nlambda * 3, 3))
  expect_equal(dim(stab$P), c(nlambda * 3, 3))
  expect_equal(dim(stab$PFER), c(nlambda * 3, 3))
  expect_equal(dim(stab$FDP), c(nlambda * 3, 3))

  # Group of outputs 2
  expect_equal(dim(stab$S_2d), c(nlambda * 3, length(pi_list), 3))

  # Group of outputs 3
  expect_equal(dim(stab$selprop), c(sum(pk), sum(pk), nlambda * 3))
  expect_equal(dim(stab$sign), c(sum(pk), sum(pk)))

  # Group of outputs 4
  default_params <- c("graphical_model", "PenalisedGraphical", "cold", "subsampling", FALSE, "MB")
  names(default_params) <- c("type", "implementation", "start", "resampling", "cpss", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$data))
  expect_equal(stab$params$lambda_other_blocks, rep(0.1, 3))

  ### Checking consistency in calibration
  expect_equal(diag(stab$Lambda[ArgmaxId(stab)[, 1], ]), Argmax(stab)[, 1])
  expect_equal(pi_list[ArgmaxId(stab)[, 2]], Argmax(stab)[, 2])
})


test_that("parallel multi-block GraphicalModel()", {
  skip_on_cran()
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- c(5, 5)
  nlambda <- 4
  K <- 10
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Data simulation
  set.seed(1)
  simul <- SimulateGraphical(n = n, pk = pk)

  stab <- GraphicalModel(
    xdata = simul$data, pk = pk,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE,
    n_cores = 2
  )

  ### Checking dimensions of the outputs
  # Group of outputs 1
  expect_equal(dim(stab$S), c(nlambda * 3, 3))
  expect_equal(dim(stab$Lambda), c(nlambda * 3, 3))
  expect_equal(dim(stab$Q), c(nlambda * 3, 3))
  expect_equal(dim(stab$Q_s), c(nlambda * 3, 3))
  expect_equal(dim(stab$P), c(nlambda * 3, 3))
  expect_equal(dim(stab$PFER), c(nlambda * 3, 3))
  expect_equal(dim(stab$FDP), c(nlambda * 3, 3))

  # Group of outputs 2
  expect_equal(dim(stab$S_2d), c(nlambda * 3, length(pi_list), 3))

  # Group of outputs 3
  expect_equal(dim(stab$selprop), c(sum(pk), sum(pk), nlambda * 3))
  expect_equal(dim(stab$sign), c(sum(pk), sum(pk)))

  # Group of outputs 4
  default_params <- c("graphical_model", "PenalisedGraphical", "warm", "subsampling", FALSE, "MB")
  names(default_params) <- c("type", "implementation", "start", "resampling", "cpss", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$data))
  expect_equal(stab$params$lambda_other_blocks, rep(0.1, 3))

  ### Checking consistency in calibration
  expect_equal(diag(stab$Lambda[ArgmaxId(stab)[, 1], ]), Argmax(stab)[, 1])
  expect_equal(pi_list[ArgmaxId(stab)[, 2]], Argmax(stab)[, 2])
})

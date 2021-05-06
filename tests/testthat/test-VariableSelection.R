test_that("outputs from VariableSelection() are of correct dimensions (gaussian)", {
  skip_on_cran()
  for (PFER_thr in c(Inf, 10)) {
    for (FDP_thr in c(Inf, 0.5)) {
      n <- 78
      pk <- 12
      nlambda <- 3
      K <- 5
      tau <- 0.55
      n_cat <- 3
      pi_list <- seq(0.6, 0.7, length.out = 15)
      simul <- SimulateRegression(n = n, pk = pk, family = "gaussian")

      if ((!is.infinite(PFER_thr)) & !is.infinite(FDP_thr)) {
        # Using only PFER_thr if both PFER_thr and FDP_thr are provided
        expect_warning(VariableSelection(
          xdata = simul$X, ydata = simul$Y,
          Lambda_cardinal = nlambda, K = K,
          pi_list = pi_list,
          tau = tau, n_cat = n_cat,
          PFER_thr = PFER_thr,
          FDP_thr = FDP_thr,
          verbose = FALSE
        ))
        stab <- suppressWarnings(VariableSelection(
          xdata = simul$X, ydata = simul$Y,
          Lambda_cardinal = nlambda, K = K,
          pi_list = pi_list,
          tau = tau, n_cat = n_cat,
          PFER_thr = PFER_thr,
          FDP_thr = FDP_thr,
          verbose = FALSE
        ))
        expect_equal(stab$params$PFER_thr, PFER_thr)
        expect_equal(stab$params$FDP_thr, Inf)
      } else {
        stab <- VariableSelection(
          xdata = simul$X, ydata = simul$Y,
          Lambda_cardinal = nlambda, K = K,
          pi_list = pi_list,
          tau = tau, n_cat = n_cat,
          PFER_thr = PFER_thr,
          FDP_thr = FDP_thr,
          verbose = FALSE
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
        expect_equal(dim(stab$selprop), c(nlambda, pk))
        expect_equal(dim(stab$Beta), c(nlambda, pk, K))

        # Group of outputs 4
        default_params <- c("glmnet", "gaussian", "subsampling", "MB")
        names(default_params) <- c("implementation", "family", "resampling", "PFER_method")
        expect_equal(unlist(stab$methods), default_params)

        # Group of outputs 5
        expect_equal(stab$params$K, K)
        expect_equal(stab$params$pi_list, pi_list)
        expect_equal(stab$params$tau, tau)
        expect_equal(stab$params$n_cat, n_cat)
        expect_equal(stab$params$pk, pk)
        expect_equal(stab$params$PFER_thr, PFER_thr)
        expect_equal(stab$params$FDP_thr, FDP_thr)
        expect_equal(stab$params$n, nrow(simul$X))
      }
    }
  }

  # Checking the error messages for incompatible glmnet families
  expect_error(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "binomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "multinomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "mgaussian"
  ))
})


test_that("outputs from VariableSelection() are of correct dimensions (binomial)", {
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- 12
  nlambda <- 3
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Binomial
  simul <- SimulateRegression(n = n, pk = pk, family = "binomial")

  stab <- VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    family = "binomial",
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE
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
  expect_equal(dim(stab$selprop), c(nlambda, pk))
  expect_equal(dim(stab$Beta), c(nlambda, pk, K))

  # Group of outputs 4
  default_params <- c("glmnet", "binomial", "subsampling", "MB")
  names(default_params) <- c("implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$X))

  # Checking the error messages for incompatible glmnet family
  # Cannot throw error for gaussian if wanted to consider 0/1 as continuous
  # Cannot throw error for multinomial as it works
  # Rely on glmnet for dealing correctly with these
  expect_error(VariableSelection(
    xdata = simul$X, ydata = simul$Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
})


test_that("outputs from VariableSelection() are of correct dimensions (mgaussian)", {
  skip_on_cran()
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- 12
  nlambda <- 3
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Binomial
  simul <- SimulateRegression(n = n, pk = pk, family = "gaussian")
  Y <- cbind(simul$Y, matrix(rnorm(nrow(simul$Y) * 2), ncol = 2))

  stab <- VariableSelection(
    xdata = simul$X, ydata = Y,
    family = "mgaussian",
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE
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
  expect_equal(dim(stab$selprop), c(nlambda, pk))
  expect_equal(dim(stab$Beta), c(nlambda, pk, K, ncol(Y)))

  # Group of outputs 4
  default_params <- c("glmnet", "mgaussian", "subsampling", "MB")
  names(default_params) <- c("implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$X))

  # Checking the error messages for incompatible glmnet families
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "gaussian"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "binomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "multinomial"
  ))
})


test_that("outputs from VariableSelection() are of correct dimensions (multinomial)", {
  skip_on_cran()
  PFER_thr <- FDP_thr <- Inf
  n <- 78
  pk <- 12
  nlambda <- 3
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Binomial
  set.seed(1)
  simul <- SimulateRegression(n = n, pk = pk, family = "binomial")
  Y <- simul$Y
  Y[Y == 0] <- sample(c(0, 2), size = sum(Y == 0), replace = TRUE)

  stab <- VariableSelection(
    xdata = simul$X, ydata = Y,
    family = "multinomial",
    Lambda_cardinal = nlambda, K = K,
    lambda.min.ratio = 0.1,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE
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
  expect_equal(dim(stab$selprop), c(nlambda, pk))
  expect_equal(dim(stab$Beta), c(nlambda, pk, K, length(unique(Y))))

  # Group of outputs 4
  default_params <- c("glmnet", "multinomial", "subsampling", "MB")
  names(default_params) <- c("implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$X))

  # Checking the error messages for incompatible glmnet families
  # Cannot throw error for gaussian if wanted to consider categories as continuous
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "binomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$X, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "mgaussian"
  ))
})

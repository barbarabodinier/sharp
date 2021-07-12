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
          xdata = simul$xdata, ydata = simul$ydata,
          Lambda_cardinal = nlambda, K = K,
          pi_list = pi_list,
          tau = tau, n_cat = n_cat,
          PFER_thr = PFER_thr,
          FDP_thr = FDP_thr,
          verbose = FALSE
        ))
        stab <- suppressWarnings(VariableSelection(
          xdata = simul$xdata, ydata = simul$ydata,
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
          xdata = simul$xdata, ydata = simul$ydata,
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
        default_params <- c("variable_selection", "PenalisedRegression", "gaussian", "subsampling", "MB")
        names(default_params) <- c("type", "implementation", "family", "resampling", "PFER_method")
        expect_equal(unlist(stab$methods), default_params)

        # Group of outputs 5
        expect_equal(stab$params$K, K)
        expect_equal(stab$params$pi_list, pi_list)
        expect_equal(stab$params$tau, tau)
        expect_equal(stab$params$n_cat, n_cat)
        expect_equal(stab$params$pk, pk)
        expect_equal(stab$params$PFER_thr, PFER_thr)
        expect_equal(stab$params$FDP_thr, FDP_thr)
        expect_equal(stab$params$n, nrow(simul$xdata))
      }
    }
  }

  # Checking the error messages for incompatible glmnet families
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "mgaussian"
  ))
})


test_that("argument penalty.factor can be used in VariableSelection()", {
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
    xdata = simul$xdata, ydata = simul$ydata,
    family = "binomial",
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE,
    penalty.factor = c(rep(1, 10), rep(0, 2))
  )
  expect_equal(ncol(stab$selprop), 10)

  # Multivariate Gaussian
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = c(6, 6), family = "gaussian")

  stab <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    family = "mgaussian",
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    PFER_thr = PFER_thr,
    FDP_thr = FDP_thr,
    verbose = FALSE,
    penalty.factor = c(rep(1, 10), rep(0, 2))
  )
  expect_equal(ncol(stab$selprop), 10)
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
    xdata = simul$xdata, ydata = simul$ydata,
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
  default_params <- c("variable_selection", "PenalisedRegression", "binomial", "subsampling", "MB")
  names(default_params) <- c("type", "implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$xdata))

  # Checking the error messages for incompatible glmnet family
  # Cannot throw error for gaussian if wanted to consider 0/1 as continuous
  # Cannot throw error for multinomial as it works
  # Rely on glmnet for dealing correctly with these
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
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

  simul <- SimulateRegression(n = n, pk = pk, family = "gaussian")
  Y <- cbind(simul$ydata, matrix(rnorm(nrow(simul$ydata) * 2), ncol = 2))

  stab <- VariableSelection(
    xdata = simul$xdata, ydata = Y,
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
  default_params <- c("variable_selection", "PenalisedRegression", "mgaussian", "subsampling", "MB")
  names(default_params) <- c("type", "implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$xdata))

  # Checking the error messages for incompatible glmnet families
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "gaussian"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "binomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
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
  n <- 200
  pk <- 12
  nlambda <- 3
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)

  # Binomial
  set.seed(1)
  simul <- SimulateRegression(n = n, pk = pk, family = "multinomial")
  Y <- simul$ydata
  Y[, 2] <- Y[, 2] * 2
  Y[, 3] <- Y[, 3] * 3
  Y <- apply(Y, 1, sum)

  stab <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
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
  default_params <- c("variable_selection", "PenalisedRegression", "multinomial", "subsampling", "MB")
  names(default_params) <- c("type", "implementation", "family", "resampling", "PFER_method")
  expect_equal(unlist(stab$methods), default_params)

  # Group of outputs 5
  expect_equal(stab$params$K, K)
  expect_equal(stab$params$pi_list, pi_list)
  expect_equal(stab$params$tau, tau)
  expect_equal(stab$params$n_cat, n_cat)
  expect_equal(stab$params$pk, pk)
  expect_equal(stab$params$PFER_thr, PFER_thr)
  expect_equal(stab$params$FDP_thr, FDP_thr)
  expect_equal(stab$params$n, nrow(simul$xdata))

  # Checking the error messages for incompatible glmnet families
  # Cannot throw error for gaussian if wanted to consider categories as continuous
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "binomial"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "cox"
  ))
  expect_error(VariableSelection(
    xdata = simul$xdata, ydata = Y,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE,
    family = "mgaussian"
  ))
})

test_that("variables with null sd in VariableSelection()", {
  n <- 78
  pk <- 12
  nlambda <- 3
  K <- 5
  tau <- 0.55
  n_cat <- 3
  pi_list <- seq(0.6, 0.7, length.out = 15)
  simul <- SimulateRegression(n = n, pk = pk, family = "gaussian")
  simul$xdata[, 1] <- rep(0, nrow(simul$xdata))

  stab <- VariableSelection(
    xdata = simul$xdata, ydata = simul$ydata,
    Lambda_cardinal = nlambda, K = K,
    pi_list = pi_list,
    tau = tau, n_cat = n_cat,
    verbose = FALSE
  )
  expect_equal(as.numeric(SelectionProportions(stab)[1]), 0)
})

test_that("cox regression in VariableSelection()", {
  # Example from glmnet
  n <- 500
  p <- 30
  nzc <- trunc(p / 10)
  beta <- rnorm(nzc)
  x <- matrix(rnorm(n * p), n, p)
  fx <- x[, seq(nzc)] %*% beta / 3
  hx <- exp(fx)
  ty <- rexp(n, hx)
  tcens <- rbinom(n = n, prob = 0.3, size = 1) # censoring indicator
  y <- cbind(time = ty, status = 1 - tcens) # y=Surv(ty,1-tcens) with library(survival)

  stab <- VariableSelection(
    xdata = x, ydata = y, K = 10, family = "cox",
    verbose = FALSE
  )
  expect_equal(as.character(stab$methods$family), "cox")
})

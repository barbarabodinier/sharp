#' Grid of penalty parameters (regression model)
#'
#' Generates a relevant grid of penalty parameter values for penalised
#' regression using the implementation in \code{\link[glmnet]{glmnet}}.
#'
#' @inheritParams VariableSelection
#' @param Lambda_cardinal number of values in the grid of parameters controlling
#'   the level of sparsity in the underlying algorithm.
#' @param check_input logical indicating if input values should be checked
#'   (recommended).
#'
#' @return A matrix of lambda values with one column and as many rows as
#'   indicated in \code{Lambda_cardinal}.
#'
#' @family lambda grid functions
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian") # simulated data
#'
#' # Lambda grid for linear regression
#' Lambda <- LambdaGridRegression(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", Lambda_cardinal = 20
#' )
#'
#' # Grid can be used in VariableSelection()
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", Lambda = Lambda
#' )
#' print(SelectedVariables(stab))
#' @export
LambdaGridRegression <- function(xdata, ydata, tau = 0.5, seed = 1,
                                 family = "gaussian",
                                 resampling = "subsampling",
                                 Lambda_cardinal = 100, check_input = TRUE,
                                 ...) {
  # Object preparation, error and warning messages
  Lambda <- NULL
  pi_list <- seq(0.6, 0.9, by = 0.01)
  K <- 100
  n_cat <- 3
  PFER_method <- "MB"
  PFER_thr <- Inf
  FDP_thr <- Inf
  verbose <- TRUE
  implementation <- PenalisedRegression
  # Checks are not re-run if coming from VariableSelection to avoid printing twice the same messages
  if (check_input) {
    # CheckInputRegression(
    #   xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
    #   K = K, tau = tau, seed = seed, n_cat = n_cat,
    #   family = family, implementation = implementation,
    #   resampling = resampling, PFER_method = PFER_method,
    #   PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    #   Lambda_cardinal = Lambda_cardinal,
    #   verbose = verbose
    # )
    # Object preparation, error and warning messages
    CheckParamRegression(
      Lambda = Lambda, pi_list = pi_list,
      K = K, tau = tau, seed = seed, n_cat = n_cat,
      family = family, implementation = implementation,
      resampling = resampling, PFER_method = PFER_method,
      PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      Lambda_cardinal = Lambda_cardinal,
      verbose = verbose
    )
    CheckDataRegression(
      xdata = xdata, ydata = ydata, family = family, verbose = verbose
    )
  }
  rm(n_cat)
  rm(Lambda)
  rm(pi_list)
  rm(K)

  # Taking one subsample/boostrap sample of the data
  withr::local_seed(1) # To keep to allow for reproducible parallelisation
  s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)

  # Applying function for variable selection to get upperbound of Lambda
  withr::local_seed(1) # To keep to allow for reproducible parallelisation
  mycv <- glmnet::glmnet(x = xdata[s, ], y = ydata[s, ], family = family, ...)
  # mycv <- do.call(glmnet::glmnet, args = list(xdata = xdata[s, ], ydata = ydata[s, ], family = family, ...))

  # Creating a grid of lambda values from min and max
  Lambda <- cbind(LambdaSequence(lmax = max(mycv$lambda), lmin = min(mycv$lambda), cardinal = Lambda_cardinal))
  Lambda <- as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda) <- paste0("s", seq(0, nrow(Lambda) - 1))

  return(Lambda)
}


#' Grid of penalty parameters (graphical model)
#'
#' Generates a relevant grid of penalty parameter values for penalised graphical
#' models.
#'
#' @inheritParams GraphicalModel
#' @param Lambda_cardinal number of values in the grid of parameters controlling
#'   the level of sparsity in the underlying algorithm.
#' @param lambda_max optional maximum value for the grid in penalty parameters.
#'   If \code{lambda_max=NULL}, the maximum value is set to the maximum
#'   covariance in absolute value. Only used if
#'   \code{implementation=PenalisedGraphical}.
#'
#' @return A matrix of lambda values with \code{length(pk)} columns and
#'   \code{Lambda_cardinal} rows.
#'
#' @family lambda grid functions
#'
#' @examples
#' \donttest{
#' # Single-block simulation
#' set.seed(1)
#' simul <- SimulateGraphical()
#'
#' # Generating grid of 10 values
#' Lambda <- LambdaGridGraphical(xdata = simul$data, Lambda_cardinal = 10)
#'
#' # Ensuring PFER < 5
#' Lambda <- LambdaGridGraphical(xdata = simul$data, Lambda_cardinal = 10, PFER_thr = 5)
#'
#' # Multi-block simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = c(10, 10))
#'
#' # Multi-block grid
#' Lambda <- LambdaGridGraphical(xdata = simul$data, pk = c(10, 10), Lambda_cardinal = 10)
#'
#' # Denser neighbouring blocks
#' Lambda <- LambdaGridGraphical(
#'   xdata = simul$data, pk = c(10, 10),
#'   Lambda_cardinal = 10, lambda_other_blocks = 0
#' )
#'
#' # Using different neighbour penalties
#' Lambda <- LambdaGridGraphical(
#'   xdata = simul$data, pk = c(10, 10),
#'   Lambda_cardinal = 10, lambda_other_blocks = c(0.1, 0, 0.1)
#' )
#' stab <- GraphicalModel(
#'   xdata = simul$data, pk = c(10, 10),
#'   Lambda = Lambda, lambda_other_blocks = c(0.1, 0, 0.1)
#' )
#' stab$Lambda
#'
#' # Visiting from empty to full graphs with max_density=1
#' Lambda <- LambdaGridGraphical(
#'   xdata = simul$data, pk = c(10, 10),
#'   Lambda_cardinal = 10, max_density = 1
#' )
#' bigblocks <- BlockMatrix(pk = c(10, 10))
#' bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
#' N_blocks <- unname(table(bigblocks_vect))
#' N_blocks # max number of edges per block
#' stab <- GraphicalModel(xdata = simul$data, pk = c(10, 10), Lambda = Lambda)
#' apply(stab$Q, 2, max, na.rm = TRUE) # max average number of edges from underlying algo
#' }
#' @export
LambdaGridGraphical <- function(xdata, pk = NULL, lambda_other_blocks = 0.1, K = 100, tau = 0.5, n_cat = 3,
                                implementation = PenalisedGraphical, start = "cold", scale = TRUE,
                                resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                                Lambda_cardinal = 50, lambda_max = NULL, lambda_path_factor = 0.001,
                                max_density = 0.5, ...) {
  # K to keep for PFER computations with PFER_method set to "SS"
  # Error and warning messages
  bigblocks <- bigblocks_vect <- blocks <- N_blocks <- nblocks <- PFER_thr_blocks <- FDP_thr_blocks <- NULL
  Lambda <- NULL
  seed <- 1 # To keep to allow for reproducible parallelisation
  pi_list <- 0.75 # only used for screening
  verbose <- TRUE
  tau_safecopy <- tau # to allow for a tau of 1 for grid definition (e.g. use with glasso without stability)

  # Need to run to define some of the objects
  CheckInputGraphical(
    xdata = xdata, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
    pi_list = pi_list, K = K, tau = 0.5, seed = seed, n_cat = n_cat,
    implementation = implementation, start = start, scale = scale,
    resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    Lambda_cardinal = Lambda_cardinal,
    lambda_max = lambda_max, lambda_path_factor = lambda_path_factor, max_density = max_density,
    verbose = verbose
  )
  tau <- tau_safecopy
  rm(Lambda)

  # Preparing lambda_dense
  ldense <- lambda_other_blocks
  p <- sum(pk)
  N <- p * (p - 1) / 2

  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in 1:ncol(xdata)) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Get upperbound of Lambda
  if (scale) {
    mycov <- stats::cor(xdata)
  } else {
    mycov <- stats::cov(xdata)
  }

  # Theoretical starting point for lambda
  if (is.null(lambda_max)) {
    diag(mycov) <- 0
    lambda_max <- max(abs(mycov))
  }
  lmin <- lambda_max
  lmin <- rep(lmin, nblocks)

  # Identifying the constraint
  if (all(is.infinite(PFER_thr))) {
    type_opt_problem <- "unconstrained"
    if (!all(is.infinite(FDP_thr))) {
      type_opt_problem <- "constrained_PFER"
      PFER_thr <- N # very loose stopping criterion for constraint on the FDP
    }
  } else {
    type_opt_problem <- "constrained_PFER"
  }

  if (type_opt_problem == "unconstrained") {
    max_q <- rep(0, nblocks)
    redo <- TRUE
    done <- rep(0, nblocks)
    while (redo) {
      lmin <- lmin * lambda_path_factor
      Lambda <- NULL
      for (b in 1:nblocks) {
        Lambda <- cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal = Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin <- Lambda[2, ]
      l <- 1
      while (l < nrow(Lambda)) {
        if (is.null(lambda_other_blocks)) {
          ldense <- lmin
        }
        tmpLambda <- Lambda[l, , drop = FALSE]
        myscreen <- SerialGraphical(
          xdata = xdata, pk = pk, Lambda = tmpLambda, lambda_other_blocks = ldense, pi_list = pi_list, K = 1,
          tau = tau, seed = seed, n_cat = n_cat,
          implementation = implementation, start = start, scale = scale,
          resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
          verbose = FALSE, ...
        ) # Only 1 iteration to get the Q

        if (l < nrow(Lambda)) {
          # Updating the smallest lambda if the density of the block is still below max_density
          for (b in 1:nblocks) {
            lmin[b] <- ifelse((myscreen$Q[b, b] < (max_density * N_blocks)[b]) & (done[b] == 0),
              yes = Lambda[l + 1, b], no = lmin[b]
            )
            done[b] <- ifelse(myscreen$Q[b, b] >= (max_density * N_blocks)[b], yes = 1, no = 0)
          }
        }

        # Increment if max_density is not yet reached
        Q_block_iteration <- NULL
        for (b in 1:nblocks) {
          Q_block_iteration <- c(Q_block_iteration, myscreen$Q[b, b])
        }

        if (any(Q_block_iteration < max_density * N_blocks)) {
          l <- l + 1
        } else {
          l <- nrow(Lambda) # stopping current while loop
          redo <- FALSE # stopping overarching while loop
        }
      }
    }
  }

  if (type_opt_problem == "constrained_PFER") {
    max_q <- rep(0, nblocks)
    redo <- TRUE
    done <- rep(0, nblocks)
    while (redo) {
      lmin <- lmin * lambda_path_factor
      Lambda <- NULL
      for (b in 1:nblocks) {
        Lambda <- cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal = Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin <- Lambda[2, ]
      l <- 1
      while (l < nrow(Lambda)) {
        if (is.null(lambda_other_blocks)) {
          ldense <- lmin
        }
        tmpLambda <- Lambda[l, , drop = FALSE]
        myscreen <- SerialGraphical(
          xdata = xdata, pk = pk, Lambda = tmpLambda, lambda_other_blocks = ldense, pi_list = pi_list, K = 1,
          tau = tau, seed = seed, n_cat = n_cat,
          resampling = resampling, scale = scale,
          implementation = implementation, start = start, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
          verbose = FALSE, ...
        ) # Only 1 iteration to get the Q

        # Compute PFER
        PFER_l <- rep(NA, nblocks)
        for (b in 1:nblocks) {
          mytmplist <- NULL
          for (j in 1:length(pi_list)) {
            pi <- pi_list[j]
            mytmplist <- c(mytmplist, PFER(q = myscreen$Q[b, b], pi = pi, N = N_blocks[b], K = K, PFER_method = PFER_method))
          }
          PFER_l[b] <- min(mytmplist)
        }

        if (l < nrow(Lambda)) {
          # Updating the smallest lambda if the PFER of the block is still below the threshold (with some margin)
          lmin <- ifelse((PFER_l <= (PFER_thr_blocks * 1.2 + 1)) & (done == 0), yes = Lambda[l + 1, ], no = lmin)
          done <- ifelse(PFER_l > (PFER_thr_blocks * 1.2 + 1), yes = 1, no = 0)
        }

        # Increment if PFER or max_density are not yet reached
        Q_block_iteration <- NULL
        for (b in 1:nblocks) {
          Q_block_iteration <- c(Q_block_iteration, myscreen$Q[b, b])
        }

        if (any(PFER_l <= (PFER_thr_blocks * 1.2 + 1)) & (any(Q_block_iteration < max_density * N_blocks))) {
          l <- l + 1
        } else {
          l <- nrow(Lambda) # stopping current while loop
          redo <- FALSE # stopping overarching while loop
        }
      }
    }
  }

  # Prepare final lambda path for each block
  Lambda <- NULL
  for (b in 1:nblocks) {
    Lambda <- cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal = Lambda_cardinal))
  }
  Lambda <- as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda) <- paste0("s", seq(0, nrow(Lambda) - 1))

  return(Lambda)
}


#' Sequence of penalty parameters
#'
#' Generates a sequence of penalty parameters from extreme values and the
#' required number of elements. The sequence is defined on the log-scale.
#'
#' @param lmax maximum value in the grid.
#' @param lmin minimum value in the grid.
#' @param cardinal number of values in the grid.
#'
#' @return A vector with values between "lmin" and "lmax" and as many values as
#'   indicated by "cardinal".
#'
#' @family lambda grid functions
#'
#' @examples
#' # Grid from extreme values
#' mygrid <- LambdaSequence(lmax = 0.7, lmin = 0.001, cardinal = 10)
#' @export
LambdaSequence <- function(lmax, lmin, cardinal = 100) {
  return(exp(seq(log(lmax), log(lmin), length.out = cardinal)))
  # return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}

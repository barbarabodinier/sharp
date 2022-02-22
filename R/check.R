#' Checking input parameters (regression model)
#'
#' Checks if input parameters are valid. For invalid parameters, this function
#' (i) stops the run and generates an error message, or (ii) sets the invalid
#' parameter to its default value and reports it in a warning message.
#'
#' @inheritParams VariableSelection
#'
#' @keywords internal
CheckInputRegression <- function(xdata, ydata = NULL, Lambda = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                                 K = 100, tau = 0.5, seed = 1, n_cat = 3,
                                 family = "gaussian", implementation = PenalisedRegression,
                                 resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                                 Lambda_cardinal = 100,
                                 verbose = TRUE) {
  # List of arguments
  myargs <- c(
    "xdata", "ydata", "Lambda", "pi_list", "K", "tau", "seed", "n_cat",
    "family",
    "PFER_method", "PFER_thr", "FDP_thr",
    "Lambda_cardinal", "verbose"
  )

  # Checking the inputs (xdata and ydata)
  xdata <- as.matrix(xdata)
  if (!is.null(ydata)) {
    if (sum(is.na(xdata)) > 0) {
      stop("Invalid input for argument 'xdata'. Missing values are not allowed in 'xdata'.")
    }
    if (sum(is.na(ydata)) > 0) {
      stop("Invalid input for argument 'ydata'. Missing values are not allowed in 'ydata'.")
    }
    if ((nrow(xdata) < 10) | (ncol(xdata) <= 1)) {
      stop("Invalid input for argument 'xdata'. Not enough data.")
    }
  }

  # Preparing xdata
  if (is.null(colnames(xdata))) {
    colnames(xdata) <- paste0("var", 1:ncol(xdata))
  }

  # Preparing ydata
  if (!is.null(ydata)) {
    if (is.vector(ydata) | is.factor(ydata)) {
      ydata <- matrix(ydata, ncol = 1)
    }
  }

  # Checking the inputs (xdata and ydata)
  if (!is.null(ydata)) {
    if (nrow(xdata) != nrow(ydata)) {
      stop("Arguments 'xdata' and 'ydata' are not compatible. They have different numbers of observations.")
    }
  }

  # Creating dummy ydata (for resampling in unsupervised models)
  if (is.null(ydata)) {
    ydata <- cbind(rep(0, nrow(xdata)))
  }

  # Naming rows of xdata and ydata
  if (is.null(rownames(xdata)) & is.null(rownames(ydata))) {
    rownames(xdata) <- paste0("obs", 1:nrow(xdata))
    rownames(ydata) <- rownames(xdata)
  } else {
    if ((is.null(rownames(xdata))) & (!is.null(rownames(ydata)))) {
      rownames(xdata) <- rownames(ydata)
    }
    if ((!is.null(rownames(xdata))) & (is.null(rownames(ydata)))) {
      rownames(ydata) <- rownames(xdata)
    }
  }

  # Re-ordering the datasets to ensure that subsamples will be the same regardless of the order of observations in the input
  ids <- sort.list(rownames(xdata))
  xdata <- xdata[ids, , drop = FALSE]
  ydata <- ydata[ids, , drop = FALSE]

  # Further checking/preparing ydata
  if ((family == "cox")) {
    if ((ncol(ydata) != 2) | (length(unique(ydata[, 2])) != 2)) {
      stop("Invalid input for argument 'ydata'. For Cox regression using glmnet, the argument 'ydata' needs to be a matrix or data frame with two columns: the time to event and binary status.")
    }
    colnames(ydata) <- c("time", "status")
    tmp <- as.factor(ydata[, 2])
    if (verbose) {
      message(paste0("Reference category: ", levels(tmp)[1]))
      message(paste0("Other category: ", levels(tmp)[2]))
    }
    ydata[, 2] <- as.numeric(tmp) - 1
    ydata <- as.matrix(ydata)
  }
  if ((family %in% c("binomial", "multinomial"))) {
    if (ncol(ydata) > 1) {
      ydata <- DummyToCategories(x = ydata, verbose = verbose)
    } else {
      ydata <- as.factor(ydata)
      if (verbose) {
        message(paste0("Reference category: ", levels(ydata)[1]))
        message(paste0("Other categorie(s): ", paste(levels(ydata)[-1], collapse = ", ")))
      }
      ydata <- as.numeric(ydata) - 1
    }
    ydata <- matrix(ydata, ncol = 1)
    rownames(ydata) <- rownames(xdata)
    ytmp <- as.numeric(table(ydata))
    if (any(ytmp == 1)) {
      stop("At least one category in 'ydata' with only one observation.")
    }
  }

  # Checking the inputs (Lambda)
  if (!is.null(Lambda)) {
    if (is.matrix(Lambda)) {
      Lambda_copy <- Lambda
      Lambda <- NULL
      for (k in 1:ncol(Lambda_copy)) {
        Lambda <- cbind(Lambda, as.numeric(Lambda_copy[, k]))
      }
    } else {
      Lambda <- as.numeric(Lambda)
      Lambda <- cbind(Lambda)
    }
    if (any(is.na(Lambda))) {
      if (all(is.na(Lambda))) {
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda <- as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
    rownames(Lambda) <- paste0("s", seq(0, nrow(Lambda) - 1))
  }

  # Checking the inputs (pi_list)
  pi_list <- sort(pi_list)
  if (n_cat == 3) {
    if (any(pi_list > 0.5) & any(pi_list < 1)) {
      if ((min(pi_list) < 0.5) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0.5) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list > 0) & any(pi_list < 1)) {
      if ((min(pi_list) < 0) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }

  # Checking the inputs (K)
  K <- as.numeric(K)
  if ((length(K) != 1) | is.na(K)) {
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K <- 100
  }

  # Checking the inputs (tau)
  tau <- as.numeric(tau)
  if ((length(tau) != 1) | is.na(tau) | (tau >= 1) | (tau <= 0)) {
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau <- 0.5
  }

  # Checking the inputs (seed)
  seed <- as.numeric(seed)
  if ((length(seed) != 1) | is.na(seed)) {
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed <- 1
  }

  # Checking the inputs (n_cat)
  n_cat <- as.numeric(n_cat)
  if ((length(n_cat) != 1) | is.na(n_cat)) {
    warning("Invalid input for argument 'n_cat'. The argument 'seed' must be set to 2 or 3. The default value (3) was used.")
    n_cat <- 3
  }

  # Checking the inputs (family)
  family <- as.character(family)
  if ((length(family) != 1) | is.na(family)) {
    stop("Invalid input for argument 'family'. The argument 'family' must be a character string.")
  }

  # Checking the inputs (implementation)
  if (!is.function(implementation)) {
    stop("Invalid input for argument 'implementation'. This argument must be a function to use for variable selection.")
  }

  # Checking the inputs (resampling)
  if ((!is.function(resampling)) & (!is.character(resampling))) {
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }

  # Checking the inputs (PFER_method)
  PFER_method <- as.character(PFER_method)
  if ((length(PFER_method) != 1) | (!PFER_method %in% c("MB", "SS"))) {
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  # Checking the inputs (PFER_method and resampling)
  if (is.character(resampling)) {
    if ((PFER_method == "SS") & (resampling == "bootstrap")) {
      warning("Arguments 'resampling' and 'PFER_method' are not compatible. With 'PFER_method' set to 'SS', the resampling is done with complementary pairs of subsamples.")
      resampling <- "subsampling"
    }
  }

  # Checking the inputs (PFER_thr)
  PFER_thr <- as.numeric(PFER_thr)
  if ((length(PFER_thr) != 1) | any(is.na(PFER_thr)) | any(PFER_thr <= 0)) {
    warning("Invalid input for argument 'PFER_thr'. The threshold in the upper-bound of the expected number of False Positives 'PFER_thr' must be a single positive number (or Inf). The default value (Inf) was used.")
    PFER_thr <- Inf
  }

  # Checking the inputs (FDP_thr)
  FDP_thr <- as.numeric(FDP_thr)
  if ((length(FDP_thr) != 1) | any(is.na(FDP_thr)) | any((!is.infinite(FDP_thr)) & (FDP_thr <= 0)) | any((!is.infinite(FDP_thr)) & (FDP_thr > 1))) {
    warning("Invalid input for argument 'FDP_thr'. The threshold in the upper-bound of the False Discovery Proportion 'FDP_thr' must be a single number between 0 and 1 (or Inf to deactivate). The default value (Inf) was used.")
    FDP_thr <- Inf
  }

  # Checking the inputs (PFER_thr and FDP_thr)
  if ((!is.infinite(PFER_thr)) & (!is.infinite(FDP_thr))) {
    warning("Arguments 'PFER_thr' and 'FDP_thr' are not compatible. Only one of these two arguments can be used (i.e. not set to Inf). Argument 'PFER_thr' was used.")
    FDP_thr <- Inf
  }

  # Checking the inputs (Lambda_cardinal)
  Lambda_cardinal <- as.numeric(Lambda_cardinal)
  if (is.null(Lambda)) {
    if ((length(Lambda_cardinal) != 1) | is.na(Lambda_cardinal) | (Lambda_cardinal < 1)) {
      warning("Invalid input for argument 'Lambda_cardinal'. The argument 'Lambda_cardinal' must be a single positive number. A value of 10 was used.")
      Lambda_cardinal <- 10
    }
  }

  # Checking the inputs (verbose)
  verbose <- as.logical(verbose)
  if ((length(verbose) != 1) | is.na(verbose)) {
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose <- TRUE
  }

  # Assigning checked values to the parent function
  for (i in 1:length(myargs)) {
    if (!is.null(get(myargs[i]))) {
      assign(myargs[i], get(myargs[i]), envir = parent.frame(n = 1))
    }
  }
}


#' Checking input parameters (graphical model)
#'
#' Checks if input parameters are valid. For invalid parameters, this function
#' (i) stops the run and generates an error message, or (ii) sets the invalid
#' parameter to its default value and reports it in a warning message.
#'
#' @inheritParams GraphicalModel
#'
#' @keywords internal
CheckInputGraphical <- function(xdata, pk = NULL, Lambda = NULL, lambda_other_blocks = 0.1,
                                pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3,
                                implementation = PenalisedGraphical, start = "cold", scale = TRUE,
                                resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                                Lambda_cardinal = 50, lambda_max = NULL, lambda_path_factor = 0.0001, max_density = 0.3,
                                verbose = TRUE) {
  # List of arguments
  myargs <- c(
    "xdata", "pk", "Lambda", "lambda_other_blocks",
    "pi_list", "K", "tau", "seed", "n_cat",
    "start", "scale",
    "PFER_method", "PFER_thr", "FDP_thr",
    "Lambda_cardinal",
    "lambda_path_factor", "max_density",
    "verbose"
  )

  # Checking the inputs (xdata)
  xdata <- as.matrix(xdata)
  if (sum(is.na(xdata)) > 0) {
    stop("Invalid input for argument 'xdata'. Missing values are not allowed in 'xdata'.")
  }
  if ((nrow(xdata) < 10) | (ncol(xdata) <= 1)) {
    stop("Invalid input for argument 'xdata'. Not enough xdata.")
  }

  # Checking the inputs (pk)
  if (!is.null(pk)) {
    pk <- as.numeric(pk)
    if (sum(pk) != ncol(xdata)) {
      stop("Invalid input for argument 'pk'. The number of variables per group 'pk' must sum to the number of columns in 'xdata'.")
    }
  } else {
    pk <- ncol(xdata)
  }

  # Checking the inputs (pi_list)
  pi_list <- sort(pi_list)
  if (n_cat == 3) {
    if (any(pi_list > 0.5) & any(pi_list < 1)) {
      if ((min(pi_list) < 0.5) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0.5) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list > 0) & any(pi_list < 1)) {
      if ((min(pi_list) < 0) | (max(pi_list) > 1)) {
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list <- pi_list[which((pi_list > 0) & (pi_list < 1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }

  # Checking the inputs (K)
  K <- as.numeric(K)
  if ((length(K) != 1) | is.na(K)) {
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K <- 100
  }

  # Checking the inputs (tau)
  tau <- as.numeric(tau)
  if ((length(tau) != 1) | is.na(tau) | (tau >= 1) | (tau <= 0)) {
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau <- 0.5
  }

  # Checking the inputs (seed)
  seed <- as.numeric(seed)
  if ((length(seed) != 1) | is.na(seed)) {
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed <- 1
  }

  # Checking the inputs (n_cat)
  n_cat <- as.numeric(n_cat)
  if ((length(n_cat) != 1) | is.na(n_cat)) {
    warning("Invalid input for argument 'n_cat'. The argument 'seed' must be set to 2 or 3. The default value (3) was used.")
    n_cat <- 3
  }

  # Checking the inputs (implementation)
  if (!is.function(implementation)) {
    stop("Invalid input for argument 'implementation'. This argument must be a function to use for graphical modelling.")
  }

  # Checking the inputs (start)
  start <- as.character(start)
  if ((length(start) != 1) | is.na(start) | (!start %in% c("cold", "warm"))) {
    warning("Invalid input for argument 'start'. The argument must be 'cold' or 'warm'. The default value (cold) was used.")
  }

  # Checking the inputs (scale)
  scale <- as.logical(scale)
  if ((length(scale) != 1) | is.na(scale)) {
    stop("Invalid input for argument 'scale'. The argument 'scale' must be logical (TRUE or FALSE).")
  }

  # Checking the inputs (resampling)
  if ((!is.function(resampling)) & (!is.character(resampling))) {
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }

  # Checking the inputs (PFER_method)
  PFER_method <- as.character(PFER_method)
  if ((length(PFER_method) != 1) | (!PFER_method %in% c("MB", "SS"))) {
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  # Checking the inputs (PFER_method and resampling)
  if (is.character(resampling)) {
    if ((PFER_method == "SS") & (resampling == "bootstrap")) {
      warning("Arguments 'resampling' and 'PFER_method' are not compatible. With 'PFER_method' set to 'SS', the resampling is done with complementary pairs of subsamples.")
      resampling <- "subsampling"
    }
  }

  # Checking the inputs (lambda_max)
  if (!is.null(lambda_max)) {
    lambda_max <- as.numeric(lambda_max)
    if ((length(lambda_max) != 1) | is.na(lambda_max) | (lambda_max <= 0)) {
      warning("Invalid input for argument 'lambda_max'. The argument 'lambda_max' must be a single positive number. The default value (NULL) was used.")
      lambda_max <- NULL
    }
  }

  # Checking the inputs (lambda_path_factor)
  lambda_path_factor <- as.numeric(lambda_path_factor)
  if ((length(lambda_path_factor) != 1) | is.na(lambda_path_factor) | (lambda_path_factor <= 0) | (lambda_path_factor >= 1)) {
    warning("Invalid input for argument 'lambda_path_factor'. The argument 'lambda_path_factor' must be a single number between 0 and 1. The default value (0.0001) was used.")
    lambda_path_factor <- 0.0001
  }

  # Checking the inputs (max_density)
  max_density <- as.numeric(max_density)
  if ((length(max_density) != 1) | is.na(max_density) | (max_density <= 0) | (max_density > 1)) {
    warning("Invalid input for argument 'max_density'. The argument 'max_density' must be a single number between 0 and 1. The default value (0.3) was used.")
    max_density <- 0.3
  }

  # Checking the inputs (Lambda_cardinal)
  Lambda_cardinal <- as.numeric(Lambda_cardinal)
  if (is.null(Lambda)) {
    if ((length(Lambda_cardinal) != 1) | is.na(Lambda_cardinal) | (Lambda_cardinal < 2)) {
      warning("Invalid input for argument 'Lambda_cardinal'. The argument 'Lambda_cardinal' must be a single positive number. A value of 10 was used.")
      Lambda_cardinal <- 10
    }
  }

  # Create matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

  # Checking the inputs (lambda_other_blocks in single-block analyses)
  if (!is.null(lambda_other_blocks)) {
    if ((length(pk) == 1)) {
      lambda_other_blocks <- NULL
    } else {
      if (length(lambda_other_blocks) == 1) {
        lambda_other_blocks <- rep(lambda_other_blocks, nblocks)
      } else {
        if (length(lambda_other_blocks) != nblocks) {
          stop(paste0(
            "Invalid input for argument 'lambda_other_blocks'. This argument must be a vector with as many entries as there are blocks in the data (i.e. ",
            nblocks, " entries in this case)."
          ))
        }
      }
    }
  }

  # Checking the inputs (verbose)
  verbose <- as.logical(verbose)
  if ((length(verbose) != 1) | is.na(verbose)) {
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose <- TRUE
  }

  # Checking the inputs (Lambda)
  if (!is.null(Lambda)) {
    if (is.matrix(Lambda)) {
      if ((ncol(Lambda) != nblocks) & (ncol(Lambda) != 1)) {
        stop(paste0("Invalid input for argument 'Lambda'. The argument 'Lambda' must be a matrix as many columns as blocks (N=", nblocks, ")."))
      }
      if (ncol(Lambda) == 1) {
        Lambda <- as.numeric(as.vector(Lambda))
      } else {
        Lambda_copy <- Lambda
        Lambda <- NULL
        for (k in 1:ncol(Lambda_copy)) {
          Lambda <- cbind(Lambda, as.numeric(Lambda_copy[, k]))
        }
      }
    } else {
      Lambda <- as.numeric(Lambda)
    }
    if (any(is.na(Lambda))) {
      if (all(is.na(Lambda))) {
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda <- as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
  }

  # Checking the inputs (PFER_thr)
  PFER_thr <- as.numeric(PFER_thr)
  if ((!length(PFER_thr) %in% c(1, nblocks)) | any(is.na(PFER_thr)) | any(PFER_thr <= 0)) {
    warning("Invalid input for argument 'PFER_thr'. The threshold in the upper-bound of the expected number of False Positives 'PFER_thr' must be a vector with positive numbers (or Inf). The default value (Inf) was used.")
    PFER_thr <- Inf
  }

  # Checking the inputs (FDP_thr)
  FDP_thr <- as.numeric(FDP_thr)
  if (length(pk) == 1) {
    if ((!length(PFER_thr) %in% c(1, nblocks)) | any(is.na(FDP_thr)) | any((!is.infinite(FDP_thr)) & (FDP_thr <= 0)) | any((!is.infinite(FDP_thr)) & (FDP_thr > 1))) {
      warning("Invalid input for argument 'FDP_thr'. The threshold in the upper-bound of the False Discovery Proportion 'FDP_thr' must be a vector with numbers between 0 and 1 (or Inf to deactivate). The default value (Inf) was used.")
      FDP_thr <- Inf
    }
  }

  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr) == 1) {
    PFER_thr_blocks <- ceiling(prop.table(N_blocks) * PFER_thr)
  } else {
    if (length(PFER_thr) == nblocks) {
      PFER_thr_blocks <- PFER_thr
    }
  }
  if (length(FDP_thr) == 1) {
    FDP_thr_blocks <- rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr) == nblocks) {
      FDP_thr_blocks <- FDP_thr
    }
  }

  # Assigning checked values to the parent function
  for (i in 1:length(myargs)) {
    assign(myargs[i], get(myargs[i]), envir = parent.frame(n = 1))
  }

  # Assigning extra objects to the parent function
  myextra <- c("bigblocks", "bigblocks_vect", "blocks", "N_blocks", "nblocks", "PFER_thr_blocks", "FDP_thr_blocks")
  for (i in 1:length(myextra)) {
    assign(myextra[i], get(myextra[i]), envir = parent.frame(n = 1))
  }
}


#' Checking that a package is installed
#'
#' Checks if a package is installed and returns an error message if not.
#'
#' @param package character string indicating the name of the package.
#'
#' @keywords internal
CheckPackageInstalled <- function(package) {
  if (!requireNamespace(package)) {
    stop(paste0("This function requires the '", package, "' package."))
  }
}

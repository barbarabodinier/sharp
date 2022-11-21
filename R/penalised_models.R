#' Penalised regression
#'
#' Runs penalised regression using implementation from
#' \code{\link[glmnet]{glmnet}}. This function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param Lambda matrix of parameters controlling the level of sparsity.
#' @param ... additional parameters passed to \code{\link[glmnet]{glmnet}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}}
#'
#' @references \insertRef{lasso}{sharp}
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mylasso <- PenalisedRegression(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = c(0.1, 0.2), family = "gaussian"
#' )
#'
#' # Using glmnet arguments
#' mylasso <- PenalisedRegression(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = c(0.1), family = "gaussian",
#'   penalty.factor = c(rep(0, 10), rep(1, 40))
#' )
#' mylasso$beta_full
#' @export
PenalisedRegression <- function(xdata, ydata, Lambda = NULL, family, ...) {
  # Checking that input data are matrices
  xdata <- as.matrix(xdata)
  ydata <- as.matrix(ydata)

  # Storing extra arguments
  extra_args <- list(...)

  # Running the regression
  if (family == "multinomial") {
    # Extracting relevant extra arguments
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "family", "type.multinomial")]

    # Running model
    mymodel <- do.call(glmnet::glmnet, args = c(list(x = xdata, y = ydata, lambda = Lambda, family = family, type.multinomial = "grouped"), tmp_extra_args))
  } else {
    # Extracting relevant extra arguments
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::glmnet)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "family")]

    # Running model
    mymodel <- do.call(glmnet::glmnet, args = c(list(x = xdata, y = ydata, lambda = Lambda, family = family), tmp_extra_args))
  }

  if (!is.infinite(mymodel$lambda[1])) {
    # Extracting and formatting the beta coefficients
    if (!family %in% c("mgaussian", "multinomial")) {
      if (length(Lambda) > 1) {
        mybeta <- suppressWarnings({
          do.call(stats::coef,
            args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
          )
        })
      } else {
        mybeta <- suppressWarnings({
          stats::coef(mymodel, s = Lambda)
        })
      }
      mybeta <- t(as.matrix(mybeta))

      # Preparing the outputs
      beta_full <- mybeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
      if ("penalty.factor" %in% names(extra_args)) {
        selected <- ifelse(mybeta[, colnames(xdata)[which(extra_args$penalty.factor != 0)], drop = FALSE] != 0, yes = 1, no = 0)
      } else {
        selected <- ifelse(mybeta[, colnames(xdata), drop = FALSE] != 0, yes = 1, no = 0)
      }
    } else {
      if (family == "mgaussian") {
        mybeta <- array(NA,
          dim = c(length(Lambda), ncol(xdata), ncol(ydata)),
          dimnames = list(1:length(Lambda), colnames(xdata), colnames(ydata))
        )
        if (length(Lambda) > 1) {
          tmpcoefs <- suppressWarnings({
            do.call(stats::coef,
              args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
            )
          })
        } else {
          tmpcoefs <- suppressWarnings({
            stats::coef(mymodel, s = Lambda)
          })
        }
        for (y_id in 1:ncol(ydata)) {
          tmpbeta <- tmpcoefs[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
        }
      }
      if (family == "multinomial") {
        y_levels <- sort(unique(ydata))
        mybeta <- array(NA,
          dim = c(length(Lambda), ncol(xdata), length(y_levels)),
          dimnames = list(1:length(Lambda), colnames(xdata), paste0("Y", y_levels))
        )
        if (length(Lambda) > 1) {
          tmpcoefs <- suppressWarnings({
            do.call(stats::coef,
              args = c(list(object = mymodel, s = Lambda, exact = TRUE, x = xdata, y = ydata), tmp_extra_args)
            )
          })
        } else {
          tmpcoefs <- suppressWarnings({
            stats::coef(mymodel, s = Lambda)
          })
        }
        for (y_id in 1:length(y_levels)) {
          tmpbeta <- tmpcoefs[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
        }
      }

      # Preparing the outputs
      if ("penalty.factor" %in% names(extra_args)) {
        selected <- ifelse(mybeta[, colnames(xdata)[which(extra_args$penalty.factor != 0)], 1, drop = FALSE] != 0, yes = 1, no = 0)
      } else {
        selected <- ifelse(mybeta[, , 1, drop = FALSE] != 0, yes = 1, no = 0)
      }
      beta_full <- mybeta
    }
  } else {
    # Returning infinite beta is the model failed
    selected <- beta_full <- Inf
  }

  return(list(selected = selected, beta_full = beta_full))
}


#' Graphical LASSO
#'
#' Runs the graphical LASSO algorithm for estimation of a Gaussian Graphical
#' Model (GGM). This function is not using stability.
#'
#' @inheritParams GraphicalModel
#' @param xdata matrix with observations as rows and variables as columns.
#' @param Lambda matrix of parameters controlling the level of sparsity.
#' @param Sequential_template logical matrix encoding the type of procedure to
#'   use for data with multiple blocks in stability selection graphical
#'   modelling. For multi-block estimation, the stability selection model is
#'   constructed as the union of block-specific stable edges estimated while the
#'   others are weakly penalised (\code{TRUE} only for the block currently being
#'   calibrated and \code{FALSE} for other blocks). Other approaches with joint
#'   calibration of the blocks are allowed (all entries are set to \code{TRUE}).
#' @param output_omega logical indicating if the estimated precision matrices
#'   should be stored and returned.
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return An array with binary and symmetric adjacency matrices along the third
#'   dimension.
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{GraphicalModel}}
#'
#' @references \insertRef{GraphicalLassoTibshirani}{sharp}
#'
#' @details The use of the procedure from Equation (4) or (5) is controlled by
#'   the argument "Sequential_template".
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical()
#'
#' # Running graphical LASSO
#' myglasso <- PenalisedGraphical(
#'   xdata = simul$data,
#'   Lambda = matrix(c(0.1, 0.2), ncol = 1)
#' )
#'
#' # Returning estimated precision matrix
#' myglasso <- PenalisedGraphical(
#'   xdata = simul$data,
#'   Lambda = matrix(c(0.1, 0.2), ncol = 1),
#'   output_omega = TRUE
#' )
#' @export
PenalisedGraphical <- function(xdata, pk = NULL, Lambda, Sequential_template = NULL,
                               scale = TRUE, start = "cold", output_omega = FALSE, ...) {
  # Checking arguments
  if (!is.matrix(Lambda)) {
    Lambda <- matrix(Lambda, ncol = 1)
  }
  if (is.null(pk)) {
    pk <- ncol(xdata)
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Create matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

  if (is.null(Sequential_template)) {
    Sequential_template <- BlockLambdaGrid(Lambda = Lambda)$Sequential_template
  }

  # Initialisation of array storing adjacency matrices
  adjacency <- array(NA, dim = c(ncol(xdata), ncol(xdata), nrow(Lambda)))
  if (output_omega) {
    omega_full <- adjacency
  }

  for (k in 1:nrow(Lambda)) {
    # Creating penalisation matrix
    if (nblocks > 1) {
      lambdamat <- bigblocks
      for (b in 1:nblocks) {
        lambdamat[bigblocks == b] <- Lambda[k, b]
      }
    } else {
      lambdamat <- Lambda[k, 1]
    }

    # Estimation of the covariance
    if (scale) {
      cov_sub <- stats::cor(xdata)
    } else {
      cov_sub <- stats::cov(xdata)
    }

    # Extracting relevant extra arguments
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glassoFast::glassoFast)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("S", "rho", "start", "w.init", "wi.init")] # avoid duplicates and args that cannot be manually set (related to warm start)

    # Estimation of the sparse inverse covariance
    if ((start == "warm") & (k != 1)) {
      if (all(which(Sequential_template[k, ]) == which(Sequential_template[k - 1, ]))) {
        # Warm start using
        g_sub <- do.call(glassoFast::glassoFast, args = c(
          list(
            S = cov_sub, rho = lambdamat,
            start = "warm", w.init = sigma, wi.init = omega
          ),
          tmp_extra_args
        ))
      } else {
        # Cold start if first iteration for the block
        g_sub <- do.call(glassoFast::glassoFast, args = c(
          list(S = cov_sub, rho = lambdamat),
          tmp_extra_args
        ))
      }
    } else {
      # Cold start
      g_sub <- do.call(glassoFast::glassoFast, args = c(
        list(S = cov_sub, rho = lambdamat),
        tmp_extra_args
      ))
    }
    omega <- g_sub$wi
    sigma <- g_sub$w

    # Creating adjacency matrix
    A <- ifelse(omega != 0, yes = 1, no = 0)
    A <- A + t(A)
    A <- ifelse(A != 0, yes = 1, no = 0)
    diag(A) <- 0

    adjacency[, , k] <- A

    if (output_omega) {
      omega_full[, , k] <- omega
    }
  }

  if (output_omega) {
    return(list(adjacency = adjacency, omega = omega_full))
  } else {
    return(list(adjacency = adjacency))
  }
}


#' Penalised Structural Equation Model
#'
#' Runs penalised Structural Equation Modelling using implementation from
#' \code{\link[regsem]{regsem}}. This function is not using stability.
#'
#' @inheritParams StructuralEquations
#' @param Lambda matrix of parameters controlling the level of sparsity.
#' @param n_convergence maximum number of attempts to convergence.
#' @param ... additional parameters passed to \code{\link[regsem]{regsem}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different regularisation parameters. Columns correspond to
#'   different parameters to estimated.} \item{beta_full}{array of model
#'   coefficients. Rows correspond to different regularisation parameters.
#'   Columns correspond to different parameters to estimated.}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}},
#'   \code{\link{LavaanMatrix}}
#'
#' @references \insertRef{RegSEM}{sharp}
#'
#' \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' \donttest{
#' # Definition of simulated effects
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' theta <- dag
#' theta[2, 4] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(theta = theta)
#'
#' # Running regularised SEM
#' mysem <- PenalisedSEM(
#'   xdata = simul$data, adjacency = dag,
#'   Lambda = LambdaSequence(lmax = 1, lmin = 0.01, cardinal = 5)
#' )
#' mysem$selected
#'
#' # Obtaining results in matrix format
#' LavaanMatrix(vect = mysem$selected[3, ], adjacency = dag)
#'
#' # Data simulation with latent variables
#' set.seed(1)
#' pk <- c(3, 2, 3)
#' simul <- SimulateStructural(
#'   n = 500,
#'   pk = pk,
#'   nu_between = 0.5,
#'   v_sign = 1,
#'   n_manifest = 3,
#'   ev_manifest = 0.95
#' )
#'
#' # Running regularised SEM
#' dag <- LayeredDAG(layers = pk, n_manifest = 3)
#' mysem <- PenalisedSEM(
#'   xdata = simul$data, adjacency = dag,
#'   Lambda = LambdaSequence(lmax = 1, lmin = 0.01, cardinal = 5)
#' )
#' mysem$selected
#' }
#' @export
PenalisedSEM <- function(xdata, adjacency, residual_covariance = NULL,
                         Lambda, n_convergence = 500, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Scaling the data (recommended by regsem)
  xdata <- scale(xdata)

  # Creating lavaan model
  model_spec <- LavaanModel(
    adjacency = adjacency,
    residual_covariance = residual_covariance,
    manifest = colnames(adjacency)[colnames(adjacency) %in% colnames(xdata)]
  )

  # Running unpenalised sem
  out_lavaan <- lavaan::sem(model = model_spec, data = xdata)

  # Initialising matrix of coefficients
  beta_full <- matrix(NA, nrow = length(Lambda), ncol = length(lavaan::coef(out_lavaan)))

  # Extracting relevant extra arguments for regsem
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = regsem::regsem)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("model", "lambda")]

  # Defining default parameters for regsem
  if (!"type" %in% names(tmp_extra_args)) {
    tmp_extra_args$type <- "lasso"
    tmp_extra_args$gradFun <- "ram"
  }
  if (!"optMethod" %in% names(tmp_extra_args)) {
    tmp_extra_args$optMethod <- "rsolnp"
  }
  if (!"pars_pen" %in% names(tmp_extra_args)) {
    tmp_extra_args$pars_pen <- "regressions"
  }

  for (k in 1:length(Lambda)) {
    # Running regularised sem
    mymodel <- suppressWarnings(try(withr::with_seed(seed = 1, code = {
      do.call(regsem::regsem, args = c(
        list(model = out_lavaan, lambda = Lambda[k]),
        tmp_extra_args
      ))
    }),
    silent = TRUE
    ))

    # Second chance to convergence
    if ((mymodel$convergence == 1) | (inherits(mymodel, "try-error"))) {
      it_conv <- 1
      while (((mymodel$convergence == 1) | (inherits(mymodel, "try-error"))) & (it_conv <= n_convergence)) {
        withr::with_seed(seed = it_conv, code = {
          sampled_lambda <- exp(stats::rnorm(
            n = 1, mean = log(Lambda[k]),
            sd = 0.5 * abs(log(Lambda[k]) - log(Lambda[ifelse(k != 1, yes = k - 1, no = k + 1)]))
          ))
        })
        mymodel <- suppressWarnings(try(withr::with_seed(seed = 1, code = {
          do.call(regsem::regsem, args = c(
            list(model = out_lavaan, lambda = sampled_lambda),
            tmp_extra_args
          ))
        }),
        silent = TRUE
        ))
        it_conv <- it_conv + 1
      }
    }

    # Storing coefficients if model convergence
    if (mymodel$convergence == 0) {
      beta_full[k, ] <- unlist(mymodel$coefficients)
    }
  }

  # Defining row and column names
  rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))
  colnames(beta_full) <- names(lavaan::coef(out_lavaan))

  # Extracting selection status
  if (length(tmp_extra_args$pars_pen) == 1) {
    if (tmp_extra_args$pars_pen == "regressions") {
      selected <- ifelse(beta_full[, !grepl("~~", colnames(beta_full))] != 0, yes = 1, no = 0)
    }
    if (tmp_extra_args$pars_pen == "loadings") { # to check for latent variables
      selected <- ifelse(beta_full[, grepl("=~", colnames(beta_full))] != 0, yes = 1, no = 0)
    }
  } else {
    selected <- ifelse(beta_full != 0, yes = 1, no = 0)
  }

  return(list(selected = selected, beta_full = beta_full))
}


#' Penalised Structural Equation Model (OpenMx)
#'
#' Runs penalised Structural Equation Modelling using implementation from
#' \code{\link[OpenMx]{OpenMx}}. This function is not using stability.
#'
#' @inheritParams StructuralEquations
#' @param Lambda matrix of parameters controlling the level of sparsity (only
#'   the minimum, maximum and length are used).
#' @param penalised optional binary matrix indicating which coefficients are
#'   regularised.
#' @param ... additional parameters passed to \code{\link[regsem]{regsem}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different regularisation parameters. Columns correspond to
#'   different parameters to estimated.} \item{beta_full}{array of model
#'   coefficients. Rows correspond to different regularisation parameters.
#'   Columns correspond to different parameters to estimated.}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}},
#'   \code{\link{OpenMxMatrix}}
#'
#' @references \insertRef{RegSEM}{sharp}
#'
#'   \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' # Definition of simulated effects
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' theta <- dag
#' theta[2, 4] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(theta = theta)
#'
#' # Running regularised SEM
#' mysem <- PenalisedOpenMx(
#'   xdata = simul$data, adjacency = dag,
#'   Lambda = seq(1, 50, by = 10)
#' )
#' mysem$selected
#'
#' # Obtaining results in matrix format
#' OpenMxMatrix(vect = mysem$selected[3, ], adjacency = dag)
#' theta # simulated is the same
#'
#' @export
PenalisedOpenMx <- function(xdata,
                            adjacency,
                            penalised = NULL,
                            residual_covariance = NULL,
                            Lambda, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Checking inputs
  if (!is.null(penalised)) {
    if (any(dim(penalised) != dim(adjacency))) {
      stop("Arguments 'adjacency' and 'penalised' are not compatible. They must be of the same dimension.")
    }
  } else {
    penalised <- adjacency
  }
  penalised <- t(penalised)

  # Scaling the data
  xdata <- scale(xdata)

  # Defining RAM matrices in OpenMx format
  ram_matrices <- OpenMxModel(
    adjacency = adjacency,
    residual_covariance = residual_covariance,
    manifest = colnames(adjacency)[colnames(adjacency) %in% colnames(xdata)]
  )

  # Defining expectation
  expectation <- OpenMx::mxExpectationRAM("A", "S", "F", "M",
    dimnames = colnames(adjacency)
  )

  # Defining fit function (maximum likelihood estimation)
  ml_function <- OpenMx::mxFitFunctionML()

  # Running unpenalised model
  model_spec <- OpenMx::mxModel(
    "Model",
    OpenMx::mxData(xdata, type = "raw"),
    ram_matrices$Amat,
    ram_matrices$Smat,
    ram_matrices$Fmat,
    ram_matrices$Mmat,
    expectation,
    ml_function
  )
  unpenalised <- OpenMx::mxRun(model_spec, silent = TRUE, suppressWarnings = TRUE)

  # Running penalised estimations
  coef_to_penalise <- as.character(stats::na.exclude(ram_matrices$Amat$labels[which(penalised == 1)]))
  mymodel <- OpenMx::mxPenaltySearch(OpenMx::mxModel(
    unpenalised,
    OpenMx::mxPenaltyLASSO(
      what = coef_to_penalise,
      name = "lasso",
      lambda = min(Lambda),
      lambda.max = max(Lambda),
      lambda.step = (max(Lambda) - min(Lambda) + 1) / length(Lambda)
    ),
    OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = 0, labels = "lambda")
  ), silent = TRUE, suppressWarnings = TRUE)$compute$steps$PS$output$detail

  # Extracting estimated coefficients
  beta_full <- as.matrix(mymodel[rev(1:nrow(mymodel)), seq(6, ncol(mymodel) - 1)])

  # Setting as missing if model did not converge
  beta_full[which(as.character(mymodel$statusCode) != "OK"), ] <- NA

  # Defining row and column names
  rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))

  # Extracting selection status
  selected <- ifelse(abs(beta_full[, coef_to_penalise, drop = FALSE]) > 1e-5, yes = 1, no = 0)

  return(list(selected = selected, beta_full = beta_full))
}


#' Writing lavaan model
#'
#' Returns model specification in \code{\link[lavaan]{lavaan}} syntax from (i)
#' the adjacency matrix of a Directed Acyclic Graph (asymmetric matrix A in
#' Reticular Action Model notation), and (ii) a binary matrix encoding nonzero
#' entries in the residual covariance matrix (symmetric matrix S in Reticular
#' Action Model notation).
#'
#' @inheritParams PenalisedSEM
#' @param manifest optional vector of manifest variable names.
#'
#' @return A character string that can be used in argument \code{model} in
#'   \code{\link[lavaan]{sem}}.
#'
#' @seealso \code{\link{PenalisedSEM}}, \code{\link{LavaanMatrix}}
#'
#' @references \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' # Definition of the model structure
#' layers <- list(
#'   c("var1", "var2", "var3"),
#'   c("var4", "var5"),
#'   c("var6", "var7", "var8")
#' )
#' dag <- LayeredDAG(layers)
#'
#' # Writing lavaan syntax
#' model_spec <- LavaanModel(adjacency = dag)
#'
#' # Checking the matrices generated by lavaan
#' mylavaan <- lavaan::sem(model = model_spec)
#' regsem::extractMatrices(mylavaan)$A
#' regsem::extractMatrices(mylavaan)$S
#'
#' # Including residual correlation
#' res_cov <- diag(ncol(dag))
#' res_cov[1, 2] <- res_cov[2, 1] <- 1
#' model_spec <- LavaanModel(
#'   adjacency = dag,
#'   residual_covariance = res_cov
#' )
#'
#' # Checking the matrices generated by lavaan
#' mylavaan <- lavaan::sem(model = model_spec)
#' regsem::extractMatrices(mylavaan)$A
#' regsem::extractMatrices(mylavaan)$S
#'
#' # Incorporating latent variables
#' dag <- LayeredDAG(layers = c(2, 1), n_manifest = 2)
#' LavaanModel(dag, manifest = paste0("x", 1:6))
#'
#' @export
LavaanModel <- function(adjacency, residual_covariance = NULL, manifest = NULL) {
  # Identifying manifest and latent variables
  if (is.null(manifest)) {
    ids_manifest <- 1:ncol(adjacency)
    ids_latent <- NULL
  } else {
    ids_manifest <- which(colnames(adjacency) %in% manifest)
    ids_latent <- which(!colnames(adjacency) %in% manifest)
  }
  # Creating residual covariance matrix structure if not provided
  if (is.null(residual_covariance)) {
    residual_covariance <- diag(ncol(adjacency))
  }

  # Checking row and column names
  if (is.null(rownames(adjacency))) {
    rownames(adjacency) <- colnames(adjacency) <- paste0("var", 1:ncol(adjacency))
  }
  rownames(residual_covariance) <- colnames(residual_covariance) <- rownames(adjacency)

  # Initialising model specification
  model_spec <- ""

  if (length(ids_latent) > 0) {
    # Listing regressions
    adjacency_latent <- adjacency[ids_latent, ids_latent]
    for (j in 1:ncol(adjacency_latent)) {
      predictors <- rownames(adjacency_latent)[which(adjacency_latent[, j] != 0)]
      if (length(predictors) > 0) {
        model_spec <- paste0(
          model_spec, colnames(adjacency_latent)[j], " ~ ",
          paste(predictors, collapse = " + "),
          " \n "
        )
      }
    }

    # Listing measurements
    for (j in ids_latent) {
      measured <- rownames(adjacency)[which(adjacency[j, ] != 0)]
      measured <- measured[measured %in% manifest]
      if (length(measured) > 0) {
        model_spec <- paste0(
          model_spec, colnames(adjacency)[j], " =~ ",
          paste(measured, collapse = " + "),
          " \n "
        )
      }
    }
  } else {
    # Listing regressions
    for (j in 1:ncol(adjacency)) {
      predictors <- rownames(adjacency)[which(adjacency[, j] != 0)]
      if (length(predictors) > 0) {
        model_spec <- paste0(
          model_spec, colnames(adjacency)[j], " ~ ",
          paste(predictors, collapse = " + "),
          " \n "
        )
      }
    }
  }

  # Listing residual correlations
  residual_covariance_upper <- residual_covariance
  residual_covariance_upper[lower.tri(residual_covariance_upper, diag = TRUE)] <- 0
  for (i in 1:(ncol(residual_covariance_upper) - 1)) {
    tmprow <- residual_covariance_upper[i, upper.tri(residual_covariance_upper)[i, ], drop = FALSE]
    correlated <- colnames(tmprow)[which(tmprow != 0)]
    independent <- colnames(tmprow)[which(tmprow == 0)]
    model_spec <- paste0(model_spec, rownames(residual_covariance_upper)[i], " ~~ ")
    if (length(independent) > 0) {
      model_spec <- paste0(
        model_spec,
        paste(paste0("0*", independent), collapse = " + ")
      )
      if (length(correlated) > 0) {
        model_spec <- paste0(
          model_spec, "+"
        )
      }
    }
    if (length(correlated) > 0) {
      model_spec <- paste0(
        model_spec,
        paste(correlated, collapse = " + ")
      )
    }
    model_spec <- paste0(model_spec, " \n ")
  }

  return(model_spec)
}


#' Matrix from lavaan outputs
#'
#' Returns a matrix from output in \code{\link[lavaan]{lavaan}} format.
#'
#' @inheritParams LavaanModel
#' @param vect vector of coefficients to assign to entries of the matrix.
#'
#' @return An asymmetric matrix.
#'
#' @seealso \code{\link{LavaanModel}}, \code{\link{PenalisedSEM}}
#'
#' @references \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' # Definition of simulated effects
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' theta <- dag
#' theta[2, 4] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(theta = theta)
#'
#' # Running regularised SEM
#' mysem <- PenalisedSEM(
#'   xdata = simul$data, adjacency = dag,
#'   Lambda = 0.3
#' )
#'
#' # Obtaining results in matrix format
#' LavaanMatrix(vect = mysem$beta_full, adjacency = dag)
#'
#' @export
LavaanMatrix <- function(vect, adjacency, residual_covariance = NULL, manifest = NULL) {
  # Creating lavaan model
  model_spec <- LavaanModel(adjacency = adjacency, residual_covariance = residual_covariance, manifest = manifest)

  # Running unpenalised sem
  out_lavaan <- lavaan::sem(model = model_spec, data = NULL)
  A <- regsem::extractMatrices(out_lavaan)$A

  # Assigning the effects to corresponding matrix entries
  for (k in 1:length(vect)) {
    A[which(A == k, arr.ind = TRUE)] <- vect[k]
  }

  # Transposing for causes as rows and consequences as columns
  A <- t(A)

  # Re-ordering as in input
  A <- A[rownames(adjacency), colnames(adjacency)]

  return(A)
}


#' Writing OpenMx model (matrix specification)
#'
#' Returns matrix specification for use in \code{\link[OpenMx]{mxModel}} from
#' (i) the adjacency matrix of a Directed Acyclic Graph (asymmetric matrix A in
#' Reticular Action Model notation), and (ii) a binary matrix encoding nonzero
#' entries in the residual covariance matrix (symmetric matrix S in Reticular
#' Action Model notation).
#'
#' @inheritParams LavaanModel
#'
#' @return A character string that can be used in argument \code{model} in
#'   \code{\link[lavaan]{sem}}.
#'
#' @seealso \code{\link{PenalisedSEM}}, \code{\link{LavaanMatrix}}
#'
#' @references \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' # Definition of simulated effects
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' theta <- dag
#' theta[2, 4] <- 0
#' theta[3, 7] <- 0
#' theta[4, 7] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(n = 500, v_between = 1, theta = theta)
#'
#' # Writing RAM matrices for mxModel
#' ram_matrices <- OpenMxModel(adjacency = dag)
#'
#' # Running unpenalised model
#' unpenalised <- OpenMx::mxRun(OpenMx::mxModel(
#'   "Model",
#'   OpenMx::mxData(simul$data, type = "raw"),
#'   ram_matrices$Amat,
#'   ram_matrices$Smat,
#'   ram_matrices$Fmat,
#'   ram_matrices$Mmat,
#'   OpenMx::mxExpectationRAM("A", "S", "F", "M", dimnames = colnames(dag)),
#'   OpenMx::mxFitFunctionML()
#' ), silent = TRUE, suppressWarnings = TRUE)
#' unpenalised$A$values
#'
#' # Incorporating latent variables
#' ram_matrices <- OpenMxModel(
#'   adjacency = dag,
#'   manifest = paste0("x", 1:7)
#' )
#' ram_matrices$Fmat$values
#'
#' # Running unpenalised model
#' unpenalised <- OpenMx::mxRun(OpenMx::mxModel(
#'   "Model",
#'   OpenMx::mxData(simul$data[, 1:7], type = "raw"),
#'   ram_matrices$Amat,
#'   ram_matrices$Smat,
#'   ram_matrices$Fmat,
#'   ram_matrices$Mmat,
#'   OpenMx::mxExpectationRAM("A", "S", "F", "M", dimnames = colnames(dag)),
#'   OpenMx::mxFitFunctionML()
#' ), silent = TRUE, suppressWarnings = TRUE)
#' unpenalised$A$values
#'
#' @export
OpenMxModel <- function(adjacency, residual_covariance = NULL, manifest = NULL) {
  # Transposing the adjacency matrix to get support of A matrix in RAM notation
  adjacency <- t(adjacency)

  # Creating residual covariance matrix structure if not provided
  if (is.null(residual_covariance)) {
    residual_covariance <- diag(ncol(adjacency))
  }

  # Creating manifest if not provided (considering all variables as measured)
  if (is.null(manifest)) {
    manifest <- colnames(adjacency)
  }
  latent <- colnames(adjacency)[!colnames(adjacency) %in% manifest]

  # Defining coefficient labels (A)
  Alabels <- matrix(paste0("a", as.vector(row(adjacency)), "_", as.vector(col(adjacency))),
    nrow = nrow(adjacency), ncol = ncol(adjacency)
  )
  Alabels <- ifelse(adjacency == 1, yes = Alabels, no = NA)

  # Defining coefficient labels (S)
  Slabels <- matrix(paste0("s", as.vector(row(residual_covariance)), "_", as.vector(col(residual_covariance))),
    nrow = nrow(residual_covariance), ncol = ncol(residual_covariance)
  )
  Slabels[lower.tri(Slabels)] <- Slabels[upper.tri(Slabels)] # symmetry
  Slabels <- ifelse(residual_covariance == 1, yes = Slabels, no = NA)

  # Creating matrix A
  adjacency_free <- ifelse(adjacency == 1, yes = TRUE, no = FALSE)
  if (length(latent) > 0) {
    for (i in 1:length(latent)) {
      id <- which(adjacency_free[, latent[i]])[1]
      adjacency_free[id, latent[i]] <- FALSE
    }
  }
  Amat <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(adjacency),
    ncol = ncol(adjacency),
    name = "A",
    free = as.vector(adjacency_free),
    values = as.vector(adjacency),
    labels = as.vector(Alabels)
  )

  # Creating matrix S
  Smat <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(residual_covariance),
    ncol = ncol(residual_covariance),
    name = "S",
    free = as.vector(ifelse(residual_covariance == 1, yes = TRUE, no = FALSE)),
    values = as.vector(residual_covariance),
    labels = as.vector(Slabels)
  )

  # Creating matrix F
  tmpmat <- diag(ncol(adjacency))
  rownames(tmpmat) <- colnames(tmpmat) <- colnames(adjacency)
  tmpmat <- tmpmat[manifest, ]
  Fmat <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(tmpmat),
    ncol = ncol(tmpmat),
    name = "F",
    free = FALSE,
    values = as.vector(tmpmat)
  )

  # Creating matrix M (all zero for scaled data)
  Mmat <- OpenMx::mxMatrix(
    type = "Full",
    nrow = 1,
    ncol = ncol(adjacency),
    name = "M",
    free = FALSE,
    values = rep(0, ncol(adjacency))
  )

  # Preparing output
  out <- list(
    Amat = Amat,
    Smat = Smat,
    Fmat = Fmat,
    Mmat = Mmat
  )

  return(out)
}


#' Matrix from OpenMx outputs
#'
#' Returns a matrix from output of \code{\link[OpenMx]{mxPenaltySearch}}.
#'
#' @inheritParams PenalisedSEM
#' @param vect vector of coefficients to assign to entries of the matrix.
#'
#' @return An asymmetric matrix.
#'
#' @seealso \code{\link{LavaanModel}}, \code{\link{PenalisedSEM}}
#'
#' @references \insertRef{lavaanBook}{sharp}
#'
#' @examples
#' # Definition of simulated effects
#' dag <- LayeredDAG(layers = c(3, 2, 3))
#' theta <- dag
#' theta[2, 4] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(theta = theta)
#'
#' # Running regularised SEM
#' mysem <- PenalisedOpenMx(
#'   xdata = simul$data, adjacency = dag,
#'   Lambda = seq(1, 50, by = 10)
#' )
#' mysem$selected
#'
#' # Obtaining results in matrix format
#' OpenMxMatrix(vect = mysem$selected[3, ], adjacency = dag)
#'
#' @export
OpenMxMatrix <- function(vect, adjacency, residual_covariance = NULL) {
  # Defining RAM matrices in OpenMx format
  ram_matrices <- OpenMxModel(adjacency = adjacency, residual_covariance = residual_covariance)

  # Assigning estimates to corresponding matrix entries
  A <- ram_matrices$Amat$values
  for (k in 1:length(vect)) {
    A[which(ram_matrices$Amat$labels == names(vect)[k], arr.ind = TRUE)] <- vect[k]
  }
  A <- t(A)

  # Defining row and column names
  rownames(A) <- rownames(adjacency)
  colnames(A) <- colnames(adjacency)

  # # Removing observed variables
  # if (!include_manifest){
  #   ids_latent=which(apply(ram_matrices$Fmat$values,2,sum)==1)
  #   A=A[ids_latent,ids_latent]
  # }

  # Defining the class of the output
  class(A) <- c("matrix", "adjacency_matrix")

  return(A)
}

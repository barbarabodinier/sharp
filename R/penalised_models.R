#' Penalised regression
#'
#' Runs penalised regression using implementation from
#' \code{\link[glmnet]{glmnet}}. This function is not using stability.
#'
#' @inheritParams VariableSelection
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
#' @examples
#'
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
  xdata=as.matrix(xdata)
  ydata=as.matrix(ydata)
  
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
        selected <- ifelse(mybeta[, colnames(xdata)[which(extra_args$penalty.factor == 1)], drop = FALSE] != 0, yes = 1, no = 0)
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
        selected <- ifelse(mybeta[, colnames(xdata)[which(extra_args$penalty.factor == 1)], 1, drop = FALSE] != 0, yes = 1, no = 0)
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


#' GGM estimation algorithm
#'
#' Runs the algorithm for estimation of a Gaussian Graphical Model (GGM)
#' specified in the argument \code{implementation} and returns the estimated
#' adjacency matrix. This function is not using stability.
#'
#' @inheritParams GraphicalModel
#' @param xdata matrix with observations as rows and variables as columns.
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
#' @details The use of the procedure from Equation (4) or (5) is controlled by
#'   the argument "Sequential_template".
#'
#' @examples
#' \dontrun{
#'
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
#' }
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

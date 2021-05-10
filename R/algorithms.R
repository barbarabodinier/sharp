#' Variable selection algorithm
#'
#' Runs the variable selection algorithm specified in the argument
#' \code{implementation} and returns matrices of model coefficients. This
#' function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mylasso <- SelectionAlgo(
#'   xdata = simul$X, ydata = simul$Y,
#'   Lambda = c(0.1, 0.2), family = "gaussian"
#' )
#'
#' # Simulation of additional outcomes
#' set.seed(2)
#' Y <- cbind(simul$Y, matrix(rnorm(nrow(simul$Y) * 2), ncol = 2))
#'
#' # Running multivariate Gaussian LASSO
#' mylasso <- SelectionAlgo(
#'   xdata = simul$X, ydata = Y,
#'   Lambda = c(0.1, 0.2), family = "mgaussian"
#' )
#' str(mylasso)
#' }
#'
#' @export
SelectionAlgo <- function(xdata, ydata, Lambda, family, implementation = PenalisedRegression, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- apply(xdata, 2, stats::sd)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }
  xdata <- scale(xdata)

  # if (implementation == "glmnet") {
  #   # Running the regression
  #   if (family == "multinomial") {
  #     mymodel <- glmnet::glmnet(x = xdata, y = ydata, lambda = Lambda, family = family, type.multinomial = "grouped", ...)
  #   } else {
  #     mymodel <- glmnet::glmnet(x = xdata, y = ydata, lambda = Lambda, family = family, ...)
  #   }
  #
  #   if (!is.infinite(mymodel$lambda[1])) {
  #     # Extracting and formatting the beta coefficients
  #     if (!family %in% c("mgaussian", "multinomial")) {
  #       mybeta <- stats::coef(mymodel)
  #       mybeta <- t(as.matrix(mybeta))
  #       mybeta <- mybeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
  #
  #       # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  #       if (any(mysd == 0)) {
  #         mybeta[, which(mysd == 0)] <- 0
  #       }
  #
  #       # Preparing the outputs
  #       selected <- ifelse(mybeta != 0, yes = 1, no = 0)
  #       beta_full <- mybeta
  #     } else {
  #       if (family == "mgaussian") {
  #         mybeta <- array(NA,
  #           dim = c(length(Lambda), ncol(xdata), ncol(ydata)),
  #           dimnames = list(paste0("s", 0:(length(Lambda) - 1)), colnames(xdata), colnames(ydata))
  #         )
  #         for (y_id in 1:ncol(ydata)) {
  #           tmpbeta <- stats::coef(mymodel)[[y_id]]
  #           tmpbeta <- t(as.matrix(tmpbeta))
  #           tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
  #           mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
  #
  #           # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  #           if (any(mysd == 0)) {
  #             mybeta[, which(mysd == 0), y_id] <- 0
  #           }
  #         }
  #       }
  #       if (family == "multinomial") {
  #         y_levels <- sort(unique(ydata))
  #         mybeta <- array(NA,
  #           dim = c(length(Lambda), ncol(xdata), length(y_levels)),
  #           dimnames = list(
  #             paste0("s", 0:(length(Lambda) - 1)), colnames(xdata),
  #             paste0("Y", y_levels)
  #           )
  #         )
  #         for (y_id in 1:length(y_levels)) {
  #           tmpbeta <- stats::coef(mymodel)[[y_id]]
  #           tmpbeta <- t(as.matrix(tmpbeta))
  #           tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
  #           mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta
  #
  #           # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  #           if (any(mysd == 0)) {
  #             mybeta[, which(mysd == 0), y_id] <- 0
  #           }
  #         }
  #       }
  #
  #       # Preparing the outputs
  #       selected <- ifelse(mybeta[, , 1, drop = FALSE] != 0, yes = 1, no = 0)
  #       beta_full <- mybeta
  #     }
  #   } else {
  #     # Returning infinite beta is the model failed
  #     selected <- beta_full <- Inf
  #   }
  # } else {
  # Applying user-defined function for variable selection
  mybeta <- do.call(implementation, args = list(xdata = xdata, ydata = ydata, Lambda = Lambda, family = family, ...))
  selected <- mybeta$selected
  beta_full <- mybeta$beta_full

  # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  if (any(mysd == 0)) {
    if (length(dim(beta_full)) == 2) {
      selected[, which(mysd == 0)] <- 0
      beta_full[, which(mysd == 0)] <- 0
    }
    if (length(dim(beta_full)) == 3) {
      selected[, which(mysd == 0)] <- 0
      beta_full[, which(mysd == 0), ] <- 0
    }
  }
  # }

  return(list(selected = selected, beta_full = beta_full))
}


#' Graph estimation algorithm
#'
#' Runs the algorithm for estimation of an undirected graph with no self-loops
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
#' myglasso <- GraphicalAlgo(xdata = simul$data, Lambda = matrix(c(0.1, 0.2), ncol = 1))
#' }
#' @export
GraphicalAlgo <- function(xdata, pk = NULL, Lambda, Sequential_template=NULL,
                          scale = TRUE, implementation = PenalisedGraphical, start = "cold", ...) {
  if (is.null(pk)) {
    pk <- ncol(xdata)
  }

  # Identifying potential variables with null standard deviation in the subsample
  mysd <- apply(xdata, 2, stats::sd)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  if (is.null(Sequential_template)){
    Sequential_template=BlockLambdaGrid(Lambda=Lambda)$Sequential_template
  }

  # # Create matrix with block indices
  # bigblocks <- BlockMatrix(pk)
  # bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  # N_blocks <- unname(table(bigblocks_vect))
  # blocks <- unique(as.vector(bigblocks_vect))
  # names(N_blocks) <- blocks
  # nblocks <- max(blocks)
  #
  # # Initialisation of array storing adjacency matrices
  # adjacency <- array(NA, dim = c(ncol(xdata), ncol(xdata), nrow(Lambda)))

  # # Going over different (sets) of penalty parameters
  # for (k in 1:nrow(Lambda)) {
  #   # Creating penalisation matrix
  #   if (nblocks > 1) {
  #     lambdamat <- bigblocks
  #     for (b in 1:nblocks) {
  #       lambdamat[bigblocks == b] <- Lambda[k, b]
  #     }
  #   } else {
  #     lambdamat <- Lambda[k, 1]
  #   }
  #
  #   if (implementation == "glassoFast") {
  #     # Estimation of the covariance
  #     if (scale) {
  #       cov_sub <- stats::cor(xdata)
  #     } else {
  #       cov_sub <- stats::cov(xdata)
  #     }
  #
  #     # Estimation of the sparse inverse covariance
  #     if ((start == "warm") & (k != 1)) {
  #       if (all(which(Sequential_template[k, ]) == which(Sequential_template[k - 1, ]))) {
  #         g_sub <- glassoFast::glassoFast(
  #           S = cov_sub, rho = lambdamat,
  #           start = "warm", w.init = sigma, wi.init = omega
  #         )
  #       } else {
  #         # Cold start if first iteration for the block
  #         g_sub <- glassoFast::glassoFast(S = cov_sub, rho = lambdamat)
  #       }
  #     } else {
  #       g_sub <- glassoFast::glassoFast(S = cov_sub, rho = lambdamat)
  #     }
  #     omega <- g_sub$wi
  #     sigma <- g_sub$w
  #
  #     # Creating adjacency matrix
  #     A <- ifelse(omega != 0, yes = 1, no = 0)
  #     A <- A + t(A)
  #     A <- ifelse(A != 0, yes = 1, no = 0)
  #     diag(A) <- 0
  #   } else {
  #     A <- do.call(implementation, args = list(xdata = xdata, Lambda = lambdamat, scale = scale, ...))
  #   }
  #
  #   # Ensuring that there is no edge for variables with always the same value (null standard deviation)
  #   if (any(mysd == 0)) {
  #     A[which(mysd == 0), ] <- 0
  #     A[, which(mysd == 0)] <- 0
  #   }
  #
  #   # Storing the estimated adjacency matrix
  #   adjacency[, , k] <- A
  # }

  adjacency=do.call(implementation, args = list(xdata=xdata, pk = pk, Lambda=Lambda, Sequential_template=Sequential_template,
                                                scale = scale, start = start, ...))

  # Ensuring that there is no edge for variables with always the same value (null standard deviation)
  for (k in 1:dim(adjacency)[3]){
    if (any(mysd == 0)) {
      adjacency[which(mysd == 0), ,k] <- 0
      adjacency[, which(mysd == 0),k] <- 0
    }
  }

  return(adjacency)
}

#' Variable selection algorithm
#'
#' Runs the variable selection algorithm specified in the argument
#' "implementation" and returns matrices of model coefficients. This function is
#' not using stability.
#'
#' @param x matrix of predictors with observations as rows and variables as
#'   columns.
#' @param y vector or matrix of outcome(s).
#' @param lambda matrix of parameters controlling the underlying feature
#'   selection algorithm specified in "implementation". With
#'   implementation="glmnet", these are penalty parameters controlling the
#'   regularised model.
#' @param family type of regression model. This argument is defined as in the
#'   \code{\link[glmnet]{glmnet}} function from the glmnet package. Possible values
#'   include "gaussian" (linear regression), "binomial" (logistic regression),
#'   "multinomial" (multinomial regression), and "cox" (survival analysis). This
#'   argument is only used with implementation="glmnet", or with functions using
#'   the family argument in the same way (see example below).
#' @param implementation name of the function to use for variable selection.
#'   With implementation="glmnet", the function \code{\link[glmnet]{glmnet}} is called.
#'   Alternatively, this argument can be a character string indicating the name
#'   of a function. The function provided must use arguments called "x", "y",
#'   "lambda" and "family" and return matrices of model coefficients (see
#'   example below).
#' @param ... additional parameters passed to the function provided in
#'   "implementation".
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mylasso <- SelectionAlgo(x = simul$X, y = simul$Y, lambda = c(0.1, 0.2), family = "gaussian")
#'
#' # Simulation of additional outcomes
#' set.seed(2)
#' Y <- cbind(simul$Y, matrix(rnorm(nrow(simul$Y) * 2), ncol = 2))
#'
#' # Running multivariate Gaussian LASSO
#' mylasso <- SelectionAlgo(x = simul$X, y = Y, lambda = c(0.1, 0.2), family = "mgaussian")
#' str(mylasso)
#' stab <- VariableSelection(xdata = simul$X, ydata = Y, family = "mgaussian")
#' @export
SelectionAlgo <- function(x, y, lambda, family, implementation = "glmnet", ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- apply(x, 2, stats::sd)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      x[, k] <- x[, k] + stats::rnorm(n = nrow(x), sd = min(mysd[mysd != 0]) / 100)
    }
  }
  x <- scale(x)

  if (implementation == "glmnet") {
    # Running the regression
    if (family == "multinomial") {
      mymodel <- glmnet::glmnet(x = x, y = y, lambda = lambda, family = family, type.multinomial = "grouped", ...)
    } else {
      mymodel <- glmnet::glmnet(x = x, y = y, lambda = lambda, family = family, ...)
    }

    if (!is.infinite(mymodel$lambda[1])) {
      # Extracting and formatting the beta coefficients
      if (!family %in% c("mgaussian", "multinomial")) {
        mybeta <- stats::coef(mymodel)
        mybeta <- t(as.matrix(mybeta))
        mybeta <- mybeta[, colnames(x)] # removing the intercept if included

        # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
        if (any(mysd == 0)) {
          mybeta[, which(mysd == 0)] <- 0
        }

        # Preparing the outputs
        selected <- ifelse(mybeta != 0, yes = 1, no = 0)
        beta_full <- mybeta
      } else {
        if (family=="mgaussian"){
          mybeta <- array(NA,
                          dim = c(length(lambda), ncol(x), ncol(y)),
                          dimnames = list(paste0("s", 0:(length(lambda) - 1)), colnames(x), colnames(y))
          )
          for (y_id in 1:ncol(y)) {
            tmpbeta <- stats::coef(mymodel)[[y_id]]
            tmpbeta <- t(as.matrix(tmpbeta))
            tmpbeta <- tmpbeta[, colnames(x), drop = FALSE] # removing the intercept if included
            mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta

            # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
            if (any(mysd == 0)) {
              mybeta[, which(mysd == 0), y_id] <- 0
            }
          }
        }
        if (family=="multinomial"){
          y_levels=sort(unique(y))
          mybeta <- array(NA,
                          dim = c(length(lambda), ncol(x), length(y_levels)),
                          dimnames = list(paste0("s", 0:(length(lambda) - 1)), colnames(x),
                                          paste0("Y", y_levels))
          )
          for (y_id in 1:length(y_levels)) {
            tmpbeta <- stats::coef(mymodel)[[y_id]]
            tmpbeta <- t(as.matrix(tmpbeta))
            tmpbeta <- tmpbeta[, colnames(x), drop = FALSE] # removing the intercept if included
            mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta

            # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
            if (any(mysd == 0)) {
              mybeta[, which(mysd == 0), y_id] <- 0
            }
          }
        }

        # Preparing the outputs
        selected <- ifelse(mybeta[, , 1, drop = FALSE] != 0, yes = 1, no = 0)
        beta_full <- mybeta
      }
    } else {
      # Returning infinite beta is the model failed
      selected <- beta_full <- Inf
    }
  } else {
    # Applying user-defined function for variable selection
    mybeta <- do.call(get(implementation), args = list(x = x, y = y, lambda = lambda, family = family, ...))
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
  }

  return(list(selected = selected, beta_full = beta_full))
}


#' Graph estimation algorithm
#'
#' Runs the algorithm for estimation of an undirected graph with no self-loops
#' specified in the argument "implementation"
#' and returns the estimated adjacency matrix.
#' This function is not using stability.
#'
#' @param x matrix with observations as rows and variables as columns.
#' @param pk vector encoding the grouping structure.
#' Only used for multi-block stability selection.
#' For this, the variables in data have to be ordered
#' by group and argument "pk" has to be a vector
#' indicating the number of variables
#' in each of the groups (see example below).
#' If pk=NULL, single-block stability selection is performed.
#' @param Lambda matrix of parameters controlling the underlying
#' feature selection algorithm specified in "implementation".
#' With implementation="glassoFast", these are penalty parameters
#' controlling the regularised model.
#' Lambda must be a matrix with as many columns as there are entries in "pk".
#' @param Sequential_template logical matrix encoding the type of procedure
#' to use for data with multiple blocks in stability selection graphical modelling.
#' For multi-block estimation, the stability selection model is
#' constructed as the union of block-specific stable edges estimated
#' while the others are weakly penalised (TRUE only for the
#' block currently being calibrated and FALSE for other blocks).
#' Other approaches with joint calibration of the blocks are allowed
#' (all entries are set to TRUE).
#' @param scale logical indicating if the correlation (if scale=TRUE)
#' or covariance (if scale=FALSE) matrix should be used as input
#' for the graphical LASSO. If implementation is not set to "glassoFast",
#' this argument must be used as input of the function provided instead.
#' @param implementation name of the function to use for graphical modelling.
#' With implementation="glassoFast", the function \code{\link[glassoFast]{glassoFast}}
#' is used for regularised estimation of a conditional independence graph.
#' Alternatively, this argument can be a character string indicating the name of a function.
#' The function provided must use arguments called "x", "lambda" and "scale"
#' and return a binary and symmetric adjacency matrix (see example below).
#' @param start character string indicating if the algorithm should be
#' initialised at the estimated (inverse) covariance with previous
#' penalty parameters (start="warm") or not (start="cold").
#' Using start="warm" can speed-up the computations.
#' Only used for implementation="glassoFast" (see argument "start"
#' in \code{\link[glassoFast]{glassoFast}}).
#' @param ... additional parameters passed to the function provided in
#' "implementation".
#'
#' @return a binary and symmetric adjacency matrix.
#'
#' @family underlying algorithm functions
#'
#' @details
#' The use of the procedure from Equation (4) or (5)
#' is controlled by the argument "Sequential_template".
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical()
#'
#' # Running graphical LASSO
#' myglasso <- GraphicalAlgo(x = simul$data, Lambda = matrix(c(0.1, 0.2), ncol = 1))
#'
#' @export
GraphicalAlgo <- function(x, pk = NULL, Lambda, Sequential_template, scale = TRUE, implementation = "glassoFast", start = "cold", ...) {
  if (is.null(pk)) {
    pk <- ncol(x)
  }

  # Identifying potential variables with null standard deviation in the subsample
  mysd <- apply(x, 2, stats::sd)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      x[, k] <- x[, k] + stats::rnorm(n = nrow(x), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Create matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

  # Initialisation of array storing adjacency matrices
  adjacency <- array(NA, dim = c(ncol(x), ncol(x), nrow(Lambda)))

  # Going over different (sets) of penalty parameters
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

    if (implementation == "glassoFast") {
      # Estimation of the covariance
      if (scale) {
        cov_sub <- stats::cor(x)
      } else {
        cov_sub <- stats::cov(x)
      }

      # Estimation of the sparse inverse covariance
      if ((start == "warm") & (k != 1)) {
        if (all(which(Sequential_template[k, ]) == which(Sequential_template[k - 1, ]))) {
          g_sub <- glassoFast::glassoFast(
            S = cov_sub, rho = lambdamat,
            start = "warm", w.init = sigma, wi.init = omega
          )
        } else {
          # Cold start if first iteration for the block
          g_sub <- glassoFast::glassoFast(S = cov_sub, rho = lambdamat)
        }
      } else {
        g_sub <- glassoFast::glassoFast(S = cov_sub, rho = lambdamat)
      }
      omega <- g_sub$wi
      sigma <- g_sub$w

      # Creating adjacency matrix
      A <- ifelse(omega != 0, yes = 1, no = 0)
      A <- A + t(A)
      A <- ifelse(A != 0, yes = 1, no = 0)
      diag(A) <- 0
    } else {
      A <- do.call(get(implementation), args = list(x = x, lambda = lambdamat, scale = scale, ...))
    }

    # Ensuring that there is no edge for variables with always the same value (null standard deviation)
    if (any(mysd == 0)) {
      A[which(mysd == 0), ] <- 0
      A[, which(mysd == 0)] <- 0
    }

    # Storing the estimated adjacency matrix
    adjacency[, , k] <- A
  }

  return(adjacency)
}

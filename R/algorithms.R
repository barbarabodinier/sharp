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
#'   Lambda = c(0.1, 0.2), family = "gaussian",
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
SelectionAlgo <- function(xdata, ydata = NULL, Lambda, family = NULL, implementation = PenalisedRegression, ...) {
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

  # Applying user-defined function for variable selection
  mybeta <- do.call(implementation, args = list(xdata = xdata, ydata = ydata, Lambda = Lambda, family = family, ...))
  selected <- mybeta$selected
  beta_full <- mybeta$beta_full

  # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  if (any(mysd == 0)) {
    selected[, which(mysd == 0)] <- 0
    if (length(dim(beta_full)) == 2) {
      beta_full[, which(mysd == 0)] <- 0
    }
    if (length(dim(beta_full)) == 3) {
      beta_full[, which(mysd == 0), ] <- 0
    }
  }

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
GraphicalAlgo <- function(xdata, pk = NULL, Lambda, Sequential_template = NULL,
                          scale = TRUE, implementation = PenalisedGraphical, start = "cold", ...) {
  if (is.null(pk)) {
    pk <- ncol(xdata)
  }

  # Identifying potential variables with null standard deviation in the subsample
  mysd <- rep(NA, ncol(xdata))
  for (j in 1:ncol(xdata)) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Setting sequence of lambda values
  if (is.null(Sequential_template)) {
    Sequential_template <- BlockLambdaGrid(Lambda = Lambda)$Sequential_template
  }

  # Computing adjacency matrices
  adjacency <- do.call(implementation, args = list(
    xdata = xdata, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
    scale = scale, start = start, ...
  ))

  # Ensuring that there is no edge for variables with always the same value (null standard deviation)
  for (k in 1:dim(adjacency)[3]) {
    if (any(mysd == 0)) {
      adjacency[which(mysd == 0), , k] <- 0
      adjacency[, which(mysd == 0), k] <- 0
    }
  }

  return(adjacency)
}

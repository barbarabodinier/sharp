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
#' # Data simulation (univariate outcome)
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mylasso <- SelectionAlgo(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = c(0.1, 0.2), family = "gaussian",
#' )
#'
#' # Data simulation (multivariate outcome)
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(15, 15, 20))
#'
#' # Running multivariate Gaussian LASSO
#' mylasso <- SelectionAlgo(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = c(0.1, 0.2), family = "mgaussian"
#' )
#' str(mylasso)
#' }
#'
#' @export
SelectionAlgo <- function(xdata, ydata = NULL,
                          Lambda, group_x = NULL,
                          family = NULL, implementation = PenalisedRegression, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in 1:ncol(xdata)) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  names(mysd) <- colnames(xdata)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Applying user-defined function for variable selection
  mybeta <- do.call(implementation, args = list(xdata = xdata, ydata = ydata, Lambda = Lambda, group_x = group_x, family = family, ...))
  selected <- mybeta$selected
  beta_full <- mybeta$beta_full

  # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
  if (!is.infinite(selected[1])) {
    if (any(mysd == 0)) {
      ids_no_sd <- intersect(names(mysd)[which(mysd == 0)], colnames(selected))
      if (length(ids_no_sd) > 0) {
        selected[, ids_no_sd] <- 0
        if (length(dim(beta_full)) == 2) {
          beta_full[, ids_no_sd] <- 0
        }
        if (length(dim(beta_full)) == 3) {
          beta_full[, ids_no_sd, ] <- 0
        }
      }
    }
  }

  return(list(selected = selected, beta_full = beta_full))
}


#' Graphical model algorithm
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
#' myglasso <- GraphicalAlgo(
#'   xdata = simul$data,
#'   Lambda = cbind(c(0.1, 0.2))
#' )
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Running graphical LASSO
#' myglasso <- GraphicalAlgo(
#'   xdata = t(simul$data),
#'   Lambda = cbind(c(2, 3)),
#'   implementation = HierarchicalClustering
#' )
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
  if ("rows" %in% names(formals(implementation))) {
    # Clustering algorithm (Lambda becomes)
    adjacency <- do.call(implementation, args = list(
      xdata = xdata, pk = pk, nc = Lambda, Sequential_template = Sequential_template,
      scale = scale, start = start, rows = FALSE, ...
    ))$comembership
  } else {
    adjacency <- do.call(implementation, args = list(
      xdata = xdata, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
      scale = scale, start = start, ...
    ))
  }

  # Ensuring that there is no edge for variables with always the same value (null standard deviation)
  for (k in 1:dim(adjacency)[3]) {
    if (any(mysd == 0)) {
      adjacency[which(mysd == 0), , k] <- 0
      adjacency[, which(mysd == 0), k] <- 0
    }
  }

  return(adjacency)
}


#' (Sparse) clustering algorithm
#'
#' Runs the (sparse) clustering algorithm specified in the argument
#' \code{implementation} and returns matrices of selected variables, variable
#' weights, and the co-membership structure. This function is not using
#' stability.
#'
#' @inheritParams Clustering
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc}
#'   is not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param Lambda vector of penalty parameters.
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{weight}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'   \item{comembership}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of 15 observations belonging to 3 groups
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100
#' )
#'
#' # Running hierarchical clustering
#' myclust <- ClusteringAlgo(
#'   xdata = simul$data, nc = 2:5,
#'   implementation = HierarchicalClustering
#' )
#'
#' # Running sparse hierarchical clustering
#' myclust <- ClusteringAlgo(
#'   xdata = simul$data,
#'   Lambda = c(1.5, 2), nc = 2:5,
#'   implementation = SparseHierarchicalClustering
#' )
#' }
#'
#' @export
ClusteringAlgo <- function(xdata,
                           Lambda = NULL, nc,
                           implementation = HierarchicalClustering, ...) {
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
  out <- do.call(implementation, args = list(xdata = xdata, nc = nc, Lambda = Lambda, ...))

  if ("weight" %in% names(out)) {
    beta_full <- out$weight

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
  } else {
    beta_full <- matrix(1, nrow = length(nc), ncol = ncol(xdata))
    rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))
    colnames(beta_full) <- colnames(xdata)
  }

  return(list(comembership = out$comembership, weight = beta_full))
}

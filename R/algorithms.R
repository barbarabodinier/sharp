#' Variable selection algorithm
#'
#' Runs the variable selection algorithm specified in the argument
#' \code{implementation}. This function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param scale logical indicating if the predictor data should be scaled.
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family wrapping functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{PenalisedRegression}},
#'   \code{\link{SparsePCA}}, \code{\link{SparsePLS}}, \code{\link{GroupPLS}},
#'   \code{\link{SparseGroupPLS}}
#'
#' @examples
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
#' simul <- SimulateRegression(pk = 50, q = 3)
#'
#' # Running multivariate Gaussian LASSO
#' mylasso <- SelectionAlgo(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = c(0.1, 0.2), family = "mgaussian"
#' )
#' str(mylasso)
#' @export
SelectionAlgo <- function(xdata, ydata = NULL,
                          Lambda, group_x = NULL, scale = TRUE,
                          family = NULL,
                          implementation = PenalisedRegression, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in seq_len(ncol(xdata))) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  names(mysd) <- colnames(xdata)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Scaling the predictor data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Applying user-defined function for variable selection
  mybeta <- do.call(implementation, args = list(xdata = xdata, ydata = ydata, Lambda = Lambda, group_x = group_x, family = family, ...))
  selected <- mybeta$selected
  beta_full <- mybeta$beta_full

  # Checking row and column names
  if (is.null(rownames(selected)) | is.null(rownames(beta_full))) {
    rownames(selected) <- rownames(beta_full) <- paste0("s", seq(0, nrow(beta_full) - 1))
  }
  if (is.null(colnames(selected)) | is.null(colnames(beta_full))) {
    colnames(selected) <- colnames(beta_full) <- paste0("coef", seq(1, ncol(beta_full)))
  }

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
#' Runs the algorithm specified in the argument \code{implementation} and
#' returns the estimated adjacency matrix. This function is not using stability.
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
#' @family wrapping functions
#' @seealso \code{\link{GraphicalModel}}, \code{\link{PenalisedGraphical}}
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
#' myglasso <- GraphicalAlgo(
#'   xdata = simul$data,
#'   Lambda = cbind(c(0.1, 0.2))
#' )
#' @export
GraphicalAlgo <- function(xdata, pk = NULL, Lambda, Sequential_template = NULL,
                          scale = TRUE, implementation = PenalisedGraphical, start = "cold", ...) {
  if (is.null(pk)) {
    pk <- ncol(xdata)
  }

  # Identifying potential variables with null standard deviation in the subsample
  mysd <- rep(NA, ncol(xdata))
  for (j in seq_len(ncol(xdata))) {
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
    # Clustering algorithm
    adjacency <- do.call(implementation, args = list(
      xdata = xdata, pk = pk, nc = Lambda, Sequential_template = Sequential_template,
      scale = scale, start = start, rows = FALSE, ...
    ))$comembership
  } else {
    adjacency <- do.call(implementation, args = list(
      xdata = xdata, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
      scale = scale, start = start, output_omega = FALSE, ...
    ))$adjacency
  }

  # Ensuring that there is no edge for variables with always the same value (null standard deviation)
  for (k in seq_len(dim(adjacency)[3])) {
    if (any(mysd == 0)) {
      adjacency[which(mysd == 0), , k] <- 0
      adjacency[, which(mysd == 0), k] <- 0
    }
  }

  return(adjacency)
}


#' (Weighted) clustering algorithm
#'
#' Runs the (weighted) clustering algorithm specified in the argument
#' \code{implementation} and returns matrices of variable weights, and the
#' co-membership structure. This function is not using stability.
#'
#' @inheritParams Clustering
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc} is
#'   not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param Lambda vector of penalty parameters.
#' @param ... additional parameters passed to the function provided in
#'   \code{implementation}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{weight}{array of model coefficients. Rows correspond to
#'   different model parameters. Columns correspond to predictors. Indices along
#'   the third dimension correspond to outcome variable(s).}
#'   \item{comembership}{array of model coefficients. Rows correspond to
#'   different model parameters. Columns correspond to predictors. Indices along
#'   the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
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
#' @export
ClusteringAlgo <- function(xdata,
                           nc = NULL,
                           eps = NULL,
                           Lambda = NULL,
                           scale = TRUE,
                           row = TRUE,
                           implementation = HierarchicalClustering, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- rep(NA, ncol(xdata))
  for (j in seq_len(ncol(xdata))) {
    mysd[j] <- stats::sd(xdata[, j])
  }
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing if clustering of columns
  if (!row) {
    xdata <- t(xdata)
  }

  # Applying user-defined function for variable selection
  out <- do.call(implementation, args = list(
    xdata = xdata,
    nc = nc,
    eps = eps,
    Lambda = Lambda,
    scale = scale,
    ...
  ))

  # Extracting the number of clusters
  if ("nc" %in% names(out)) {
    nc <- out$nc
  } else {
    if (is.null(Lambda)) {
      nc <- nc
    } else {
      nc <- cbind(rep(nc, nrow(Lambda)))
    }
  }

  # Extracting the noise status
  if ("noise" %in% names(out)) {
    noise <- out$noise
  } else {
    noise <- rep(NA, nrow(xdata))
  }

  if ("weight" %in% names(out)) {
    beta_full <- out$weight

    # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
    if (any(mysd == 0)) {
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

  return(list(comembership = out$comembership, weight = beta_full, nc = nc, noise = noise))
}

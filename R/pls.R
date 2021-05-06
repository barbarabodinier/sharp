#' Sparse Partial Least Squares
#'
#' Runs a sparse Partial Least Squares model using implementation from
#' \code{\link[sgPLS]{sgPLS-package}}. This function is not using stability.
#'
#' @inheritParams SelectionAlgo
#' @param lambda matrix of parameters controlling the number of selected
#'   predictors at current component, as defined by \code{ncomp}.
#' @param family type of PLS model. If \code{family="gaussian"}, a sparse PLS
#'   model as defined in \code{\link[sgPLS]{sPLS}} is run (for continuous
#'   outcomes). If \code{family="binomial"}, a PLS-DA model as defined in
#'   \code{\link[sgPLS]{sPLSda}} is run (for categorical outcomes).
#' @param ncomp number of components.
#' @param keepX_previous number of selected predictors in previous components.
#'   Only used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{lambda}.
#' @param keepY number of selected outcome variables. This argument is defined
#'   as in \code{\link[sgPLS]{sPLS}}. Only used if \code{family="gaussian"}.
#' @param ... additional arguments to be passed to \code{\link[sgPLS]{sPLS}} or
#'   \code{\link[sgPLS]{sPLSda}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors (starting
#'   with "X") or outcomes (starting with "Y") variables for different
#'   components (denoted by "PC").}
#'
#' @family PLS functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}}
#'
#' @examples
#' ## Sparse PLS
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' ydata <- cbind(simul$Y, matrix(rnorm(100 * 3), ncol = 3))
#' colnames(ydata) <- paste0("outcome", 1:4)
#' x <- simul$X
#' y <- ydata
#'
#' # Running sPLS with 2 X-variables and 1 Y-variable
#' mypls <- SparsePLS(x = x, y = y, lambda = 2, family = "gaussian", keepY = 1)
#'
#' ## Sparse PLS-DA
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#'
#' # Running sPLS-DA with 2 X-variables and 1 Y-variable
#' mypls <- SparsePLS(x = simul$X, y = simul$Y, lambda = 2, family = "binomial")
#' @export
SparsePLS <- function(x, y, lambda, family = "gaussian", ncomp = 1, keepX_previous = NULL, keepY = NULL, ...) {
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Re-formatting y
  if (is.vector(y)) {
    y <- cbind(y)
  }

  # All Y variables are selected by default
  if (is.null(keepY)) {
    keepY <- rep(ncol(y), ncomp)
  }

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(x), ncomp - 1)
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(lambda), ncol = ncol(x))
  rownames(beta) <- paste0("s", 0:(length(lambda) - 1))
  colnames(beta) <- colnames(x)

  # Initialising the full set of loadings coefficients
  if (family == "gaussian") {
    beta_full <- matrix(NA, nrow = length(lambda), ncol = (ncol(x) + ncol(y)) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(x), "_PC"), rep(1:ncomp, each = ncol(x))),
      paste0(paste0("Y_", colnames(y), "_PC"), rep(1:ncomp, each = ncol(y)))
    )
  } else {
    ncat <- length(unique(y))
    beta_full <- matrix(NA, nrow = length(lambda), ncol = (ncol(x) + ncat) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(x), "_PC"), rep(1:ncomp, each = ncol(x))),
      paste0(paste0("Y_", sort(unique(y)), "_PC"), rep(1:ncomp, each = ncat))
    )
  }

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, lambda[k])

    if (family == "gaussian") {
      mymodel <- sgPLS::sPLS(X = x, Y = y, ncomp = ncomp, keepX = nvarx, keepY = keepY, ...)
    } else {
      mymodel <- sgPLS::sPLSda(X = x, Y = as.vector(y), ncomp = ncomp, keepX = nvarx, ...) # no sparsity in Y
    }

    # Extracting X and Y loadings
    Xloadings <- mymodel$loadings$X
    Yloadings <- mymodel$loadings$Y

    # Making sure that the loadings of the univariate Y are positives
    if (ncol(y) == 1) {
      for (j in 1:ncol(mymodel$loadings$Y)) {
        if (sign(mymodel$loadings$Y[1, j])) {
          Yloadings[1, j] <- -Yloadings[1, j]
          Xloadings[, j] <- -Xloadings[, j]
        }
      }
    }

    # Exporting the set of loadings to look at for selection with current component
    beta[k, ] <- Xloadings[, ncomp, drop = FALSE]

    # Exporting all loadings coefficients
    beta_full[k, ] <- c(as.vector(Xloadings), as.vector(Yloadings))
  }

  beta <- ifelse(beta != 0, yes = 1, no = 0)

  return(list(selected = beta, beta_full = beta_full))
}


#' Sparse group Partial Least Squares
#'
#' Runs a sparse group Partial Least Squares model using implementation from
#' \code{\link[sgPLS]{sgPLS-package}}. This function is not using stability.
#'
#' @inheritParams SparsePLS
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group.
#' @param group_y optional vector encoding the grouping structure among
#'   outcomes. This argument indicates the number of variables in each group.
#' @param lambda matrix of parameters controlling the number of selected groups
#'   at current component, as defined by \code{ncomp}.
#' @param keepX_previous number of selected groups in previous components. Only
#'   used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sgPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{lambda}.
#' @param keepY number of selected groups of outcome variables. This argument is
#'   defined as in \code{\link[sgPLS]{sgPLS}}. Only used if
#'   \code{family="gaussian"}.
#' @param ... additional arguments to be passed to \code{\link[sgPLS]{sgPLS}} or
#'   \code{\link[sgPLS]{sgPLSda}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors (starting
#'   with "X") or outcomes (starting with "Y") variables for different
#'   components (denoted by "PC").}
#'
#' @family PLS functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}}
#'
#' @examples
#' ## Sparse group PLS
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' ydata <- cbind(simul$Y, matrix(rnorm(100 * 3), ncol = 3))
#' colnames(ydata) <- paste0("outcome", 1:4)
#' x <- simul$X
#' y <- ydata
#'
#' # Running sgPLS with 1 group and sparsity of 0.5
#' mypls <- SparseGroupPLS(
#'   x = x, y = y, lambda = 1, family = "gaussian",
#'   group_x = c(20, 15, 25), alpha.x = 0.5
#' )
#'
#' # Running sgPLS with groups on outcomes
#' mypls <- SparseGroupPLS(
#'   x = x, y = y, lambda = 1, family = "gaussian",
#'   group_x = c(20, 15, 25), alpha.x = 0.5,
#'   group_y = c(2, 2), keepY = 1, alpha.y = 0.9
#' )
#'
#' ## Sparse group PLS-DA
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "binomial")
#'
#' # Running sgPLS-DA with 1 group and sparsity of 0.9
#' mypls <- SparseGroupPLS(
#'   x = simul$X, y = simul$Y, lambda = 1, family = "binomial",
#'   group_x = c(20, 15, 25), alpha.x = 0.9
#' )
#' @export
SparseGroupPLS <- function(x, y, group_x, group_y = NULL,
                           lambda, family = "gaussian", ncomp = 1,
                           keepX_previous = NULL, keepY = NULL, ...) {
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Re-formatting y
  if (is.vector(y)) {
    y <- cbind(y)
  }

  # All Y variables are selected by default
  if (is.null(keepY)) {
    keepY <- rep(ncol(y), ncomp)
  }

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(x), ncomp - 1)
  }

  # Re-formatting the grouping
  ind.block.x <- cumsum(group_x)[-length(group_x)]
  if (!is.null(group_y)) {
    ind.block.y <- cumsum(group_y)[-length(group_y)]
  } else {
    ind.block.y <- NULL
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(lambda), ncol = ncol(x))
  rownames(beta) <- paste0("s", 0:(length(lambda) - 1))
  colnames(beta) <- colnames(x)

  # Initialising the full set of loadings coefficients
  if (family == "gaussian") {
    beta_full <- matrix(NA, nrow = length(lambda), ncol = (ncol(x) + ncol(y)) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(x), "_PC"), rep(1:ncomp, each = ncol(x))),
      paste0(paste0("Y_", colnames(y), "_PC"), rep(1:ncomp, each = ncol(y)))
    )
  } else {
    ncat <- length(unique(y))
    beta_full <- matrix(NA, nrow = length(lambda), ncol = (ncol(x) + ncat) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(x), "_PC"), rep(1:ncomp, each = ncol(x))),
      paste0(paste0("Y_", sort(unique(y)), "_PC"), rep(1:ncomp, each = ncat))
    )
  }

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, lambda[k])

    if (family == "binomial") {
      mymodel <- sgPLS::sgPLSda(
        X = x, Y = as.vector(y),
        ind.block.x = ind.block.x,
        ncomp = ncomp, keepX = nvarx, ...
      ) # no sparsity in Y
    } else {
      mymodel <- sgPLS::sgPLS(
        X = x, Y = y,
        ind.block.x = ind.block.x, ind.block.y = ind.block.y,
        ncomp = ncomp, keepX = nvarx, keepY = keepY, ...
      )
    }

    # Extracting X and Y loadings
    Xloadings <- mymodel$loadings$X
    Yloadings <- mymodel$loadings$Y

    # Making sure that the loadings of the univariate Y are positives
    if (ncol(y) == 1) {
      for (j in 1:ncol(mymodel$loadings$Y)) {
        if (sign(mymodel$loadings$Y[1, j])) {
          Yloadings[1, j] <- -Yloadings[1, j]
          Xloadings[, j] <- -Xloadings[, j]
        }
      }
    }

    # Exporting the set of loadings to look at for selection with current component
    beta[k, ] <- Xloadings[, ncomp, drop = FALSE]

    # Exporting all loadings coefficients
    beta_full[k, ] <- c(as.vector(Xloadings), as.vector(Yloadings))
  }

  beta <- ifelse(beta != 0, yes = 1, no = 0)

  return(list(selected = beta, beta_full = beta_full))
}

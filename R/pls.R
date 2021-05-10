#' Sparse Partial Least Squares
#'
#' Runs a sparse Partial Least Squares model using implementation from
#' \code{\link[sgPLS]{sgPLS-package}}. This function is not using stability.
#'
#' @inheritParams SelectionAlgo
#' @param Lambda matrix of parameters controlling the number of selected
#'   predictors at current component, as defined by \code{ncomp}.
#' @param family type of PLS model. If \code{family="gaussian"}, a sparse PLS
#'   model as defined in \code{\link[sgPLS]{sPLS}} is run (for continuous
#'   outcomes). If \code{family="binomial"}, a PLS-DA model as defined in
#'   \code{\link[sgPLS]{sPLSda}} is run (for categorical outcomes).
#' @param ncomp number of components.
#' @param keepX_previous number of selected predictors in previous components.
#'   Only used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{Lambda}.
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
#' mypls <- SparsePLS(xdata = x, ydata = y, Lambda = 2, family = "gaussian", keepY = 1)
#'
#' ## Sparse PLS-DA
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#'
#' # Running sPLS-DA with 2 X-variables and 1 Y-variable
#' mypls <- SparsePLS(xdata = simul$X, ydata = simul$Y, Lambda = 2, family = "binomial")
#' @export
SparsePLS <- function(xdata, ydata, Lambda, family = "gaussian", ncomp = 1, keepX_previous = NULL, keepY = NULL, ...) {
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Re-formatting y
  if (is.vector(ydata)) {
    ydata <- cbind(ydata)
  }

  # All Y variables are selected by default
  if (is.null(keepY)) {
    keepY <- rep(ncol(ydata), ncomp)
  }

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(xdata), ncomp - 1)
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(Lambda), ncol = ncol(xdata))
  rownames(beta) <- paste0("s", 0:(length(Lambda) - 1))
  colnames(beta) <- colnames(xdata)

  # Initialising the full set of loadings coefficients
  if (family == "gaussian") {
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncol(ydata)) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", colnames(ydata), "_PC"), rep(1:ncomp, each = ncol(ydata)))
    )
  } else {
    ncat <- length(unique(ydata))
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncat) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", sort(unique(ydata)), "_PC"), rep(1:ncomp, each = ncat))
    )
  }

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(Lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, Lambda[k])

    if (family == "gaussian") {
      mymodel <- sgPLS::sPLS(X = xdata, Y = ydata, ncomp = ncomp, keepX = nvarx, keepY = keepY, ...)
    } else {
      mymodel <- sgPLS::sPLSda(X = xdata, Y = as.vector(ydata), ncomp = ncomp, keepX = nvarx, ...) # no sparsity in Y
    }

    # Extracting X and Y loadings
    Xloadings <- mymodel$loadings$X
    Yloadings <- mymodel$loadings$Y

    # Making sure that the loadings of the univariate Y are positives
    if (ncol(ydata) == 1) {
      for (j in 1:ncol(mymodel$loadings$Y)) {
        if (sign(mymodel$loadings$Y[1, j]) == -1) {
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
#' @param family type of PLS model. If \code{family="gaussian"}, a sparse group
#'   PLS model as defined in \code{\link[sgPLS]{sgPLS}} is run (for continuous
#'   outcomes). If \code{family="binomial"}, a PLS-DA model as defined in
#'   \code{\link[sgPLS]{sgPLSda}} is run (for categorical outcomes).
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group.
#' @param group_y optional vector encoding the grouping structure among
#'   outcomes. This argument indicates the number of variables in each group.
#' @param alpha.x vector of parameters controlling the level of sparsity within
#'   groups of predictors.
#' @param alpha.y optional vector of parameters controlling the level of
#'   sparsity within groups of outcomes. Only used if \code{family="gaussian"}.
#' @param Lambda matrix of parameters controlling the number of selected groups
#'   at current component, as defined by \code{ncomp}.
#' @param keepX_previous number of selected groups in previous components. Only
#'   used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sgPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{Lambda}.
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
#'   xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'   group_x = c(10, 15, 25), alpha.x = 0.5
#' )
#'
#' # Running sgPLS with groups on outcomes
#' mypls <- SparseGroupPLS(
#'   xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'   group_x = c(10, 15, 25), alpha.x = 0.5,
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
#'   xdata = simul$X, ydata = simul$Y, Lambda = 1, family = "binomial",
#'   group_x = c(10, 15, 25), alpha.x = 0.9
#' )
#' @export
SparseGroupPLS <- function(xdata, ydata, family = "gaussian", group_x, group_y = NULL,
                           Lambda, alpha.x, alpha.y = NULL,
                           keepX_previous = NULL, keepY = NULL, ncomp = 1, ...) {
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Re-formatting y
  if (is.vector(ydata)) {
    ydata <- cbind(ydata)
  }

  # All Y variables are selected by default
  if (is.null(keepY)) {
    keepY <- rep(ncol(ydata), ncomp)
  }

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(xdata), ncomp - 1)
  }

  # Re-formatting the grouping
  ind.block.x <- cumsum(group_x)[-length(group_x)]
  if (!is.null(group_y)) {
    ind.block.y <- cumsum(group_y)[-length(group_y)]
  } else {
    ind.block.y <- NULL
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(Lambda), ncol = ncol(xdata))
  rownames(beta) <- paste0("s", 0:(length(Lambda) - 1))
  colnames(beta) <- colnames(xdata)

  # Initialising the full set of loadings coefficients
  if (family == "gaussian") {
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncol(ydata)) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", colnames(ydata), "_PC"), rep(1:ncomp, each = ncol(ydata)))
    )
  } else {
    ncat <- length(unique(ydata))
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncat) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", sort(unique(ydata)), "_PC"), rep(1:ncomp, each = ncat))
    )
  }

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(Lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, Lambda[k])

    if (family == "binomial") {
      mymodel <- sgPLS::sgPLSda(
        X = xdata, Y = as.vector(ydata),
        ind.block.x = ind.block.x,
        ncomp = ncomp, keepX = nvarx, alpha.x = alpha.x, ...
      ) # no sparsity in Y
    } else {
      mymodel <- sgPLS::sgPLS(
        X = xdata, Y = ydata,
        ind.block.x = ind.block.x, ind.block.y = ind.block.y,
        alpha.x = alpha.x, alpha.y = alpha.y,
        ncomp = ncomp, keepX = nvarx, keepY = keepY, ...
      )
    }

    # Extracting X and Y loadings
    Xloadings <- mymodel$loadings$X
    Yloadings <- mymodel$loadings$Y

    # Making sure that the loadings of the univariate Y are positives
    if (ncol(ydata) == 1) {
      for (j in 1:ncol(mymodel$loadings$Y)) {
        if (sign(mymodel$loadings$Y[1, j]) == -1) {
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


#' Group Partial Least Squares
#'
#' Runs a group Partial Least Squares model using implementation from
#' \code{\link[sgPLS]{sgPLS-package}}. This function is not using stability.
#'
#' @inheritParams SparsePLS
#' @param family type of PLS model. If \code{family="gaussian"}, a group PLS
#'   model as defined in \code{\link[sgPLS]{gPLS}} is run (for continuous
#'   outcomes). If \code{family="binomial"}, a PLS-DA model as defined in
#'   \code{\link[sgPLS]{gPLSda}} is run (for categorical outcomes).
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group.
#' @param group_y optional vector encoding the grouping structure among
#'   outcomes. This argument indicates the number of variables in each group.
#' @param Lambda matrix of parameters controlling the number of selected groups
#'   at current component, as defined by \code{ncomp}.
#' @param keepX_previous number of selected groups in previous components. Only
#'   used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sgPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{Lambda}.
#' @param keepY number of selected groups of outcome variables. This argument is
#'   defined as in \code{\link[sgPLS]{sgPLS}}. Only used if
#'   \code{family="gaussian"}.
#' @param ... additional arguments to be passed to \code{\link[sgPLS]{gPLS}} or
#'   \code{\link[sgPLS]{gPLSda}}.
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
#' ## Group PLS
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' ydata <- cbind(simul$Y, matrix(rnorm(100 * 3), ncol = 3))
#' colnames(ydata) <- paste0("outcome", 1:4)
#' x <- simul$X
#' y <- ydata
#'
#' # Running gPLS with 1 group and sparsity of 0.5
#' mypls <- GroupPLS(
#'   xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'   group_x = c(10, 15, 25),
#' )
#'
#' # Running gPLS with groups on outcomes
#' mypls <- GroupPLS(
#'   xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'   group_x = c(10, 15, 25),
#'   group_y = c(2, 2), keepY = 1
#' )
#'
#' ## Group PLS-DA
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "binomial")
#'
#' # Running sgPLS-DA with 1 group and sparsity of 0.9
#' mypls <- GroupPLS(
#'   xdata = simul$X, ydata = simul$Y, Lambda = 1, family = "binomial",
#'   group_x = c(10, 15, 25)
#' )
#' @export
GroupPLS <- function(xdata, ydata, family = "gaussian", group_x, group_y = NULL,
                     Lambda, keepX_previous = NULL, keepY = NULL, ncomp = 1, ...) {
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Re-formatting y
  if (is.vector(ydata)) {
    ydata <- cbind(ydata)
  }

  # All Y variables are selected by default
  if (is.null(keepY)) {
    keepY <- rep(ncol(ydata), ncomp)
  }

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(xdata), ncomp - 1)
  }

  # Re-formatting the grouping
  ind.block.x <- cumsum(group_x)[-length(group_x)]
  if (!is.null(group_y)) {
    ind.block.y <- cumsum(group_y)[-length(group_y)]
  } else {
    ind.block.y <- NULL
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(Lambda), ncol = ncol(xdata))
  rownames(beta) <- paste0("s", 0:(length(Lambda) - 1))
  colnames(beta) <- colnames(xdata)

  # Initialising the full set of loadings coefficients
  if (family == "gaussian") {
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncol(ydata)) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", colnames(ydata), "_PC"), rep(1:ncomp, each = ncol(ydata)))
    )
  } else {
    ncat <- length(unique(ydata))
    beta_full <- matrix(NA, nrow = length(Lambda), ncol = (ncol(xdata) + ncat) * ncomp)
    rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
    colnames(beta_full) <- c(
      paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))),
      paste0(paste0("Y_", sort(unique(ydata)), "_PC"), rep(1:ncomp, each = ncat))
    )
  }

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(Lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, Lambda[k])

    if (family == "binomial") {
      mymodel <- sgPLS::gPLSda(
        X = xdata, Y = as.vector(ydata),
        ind.block.x = ind.block.x,
        ncomp = ncomp, keepX = nvarx, ...
      ) # no sparsity in Y
    } else {
      mymodel <- sgPLS::gPLS(
        X = xdata, Y = ydata,
        ind.block.x = ind.block.x, ind.block.y = ind.block.y,
        ncomp = ncomp, keepX = nvarx, keepY = keepY, ...
      )
    }

    # Extracting X and Y loadings
    Xloadings <- mymodel$loadings$X
    Yloadings <- mymodel$loadings$Y

    # Making sure that the loadings of the univariate Y are positives
    if (ncol(ydata) == 1) {
      for (j in 1:ncol(mymodel$loadings$Y)) {
        if (sign(mymodel$loadings$Y[1, j]) == -1) {
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

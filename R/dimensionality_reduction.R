#' Partial Least Squares 'a la carte'
#'
#' Runs a Partial Least Squares (PLS) model in regression mode using algorithm
#' implemented in \code{\link[mixOmics]{pls}}. This function allows for the
#' construction of components based on different sets of predictor and/or
#' outcome variables. This function is not using stability.
#'
#' @inheritParams SelectionAlgo
#' @param selectedX binary matrix of size \code{(ncol(xdata) * ncomp)}. The
#'   binary entries indicate which predictors (in rows) contribute to the
#'   definition of each component (in columns). If \code{selectedX=NULL}, all
#'   predictors are selected for all components.
#' @param selectedY binary matrix of size \code{(ncol(ydata) * ncomp)}. The
#'   binary entries indicate which outcomes (in rows) contribute to the
#'   definition of each component (in columns). If \code{selectedY=NULL}, all
#'   outcomes are selected for all components.
#' @param family type of PLS model. Only \code{family="gaussian"} is supported.
#'   This corresponds to a PLS model as defined in \code{\link[mixOmics]{pls}}
#'   (for continuous outcomes).
#' @param ncomp number of components.
#' @param scale logical indicating if the data should be scaled (i.e.
#'   transformed so that all variables have a standard deviation of one).
#'
#' @details All matrices are defined as in (Wold et al. 2001). The weight matrix
#'   \code{Wmat} is the equivalent of \code{loadings$X} in
#'   \code{\link[mixOmics]{pls}}. The loadings matrix \code{Pmat} is the
#'   equivalent of \code{mat.c} in \code{\link[mixOmics]{pls}}. The score
#'   matrices \code{Tmat} and \code{Qmat} are the equivalent of
#'   \code{variates$X} and \code{variates$Y} in \code{\link[mixOmics]{pls}}.
#'
#' @return A list with: \item{Wmat}{matrix of X-weights.} \item{Wstar}{matrix of
#'   transformed X-weights.} \item{Pmat}{matrix of X-loadings.}
#'   \item{Cmat}{matrix of Y-weights.} \item{Tmat}{matrix of X-scores.}
#'   \item{Umat}{matrix of Y-scores.} \item{Qmat}{matrix needed for
#'   predictions.} \item{Rmat}{matrix needed for predictions.}
#'   \item{meansX}{vector used for centering of predictors, needed for
#'   predictions.} \item{sigmaX}{vector used for scaling of predictors, needed
#'   for predictions.} \item{meansY}{vector used for centering of outcomes,
#'   needed for predictions.} \item{sigmaY}{vector used for scaling of outcomes,
#'   needed for predictions.} \item{methods}{a list with \code{family} and
#'   \code{scale} values used for the run.} \item{params}{a list with
#'   \code{selectedX} and \code{selectedY} values used for the run.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @references \insertRef{PLS}{sharp}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   oldpar <- par(no.readonly = TRUE)
#'
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 200, pk = 15, q = 3, family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'
#'   # PLS
#'   mypls <- PLS(xdata = x, ydata = y, ncomp = 3)
#'
#'   if (requireNamespace("sgPLS", quietly = TRUE)) {
#'     # Sparse PLS to identify relevant variables
#'     stab <- BiSelection(
#'       xdata = x, ydata = y,
#'       family = "gaussian", ncomp = 3,
#'       LambdaX = 1:(ncol(x) - 1),
#'       LambdaY = 1:(ncol(y) - 1),
#'       implementation = SparsePLS,
#'       n_cat = 2
#'     )
#'     plot(stab)
#'
#'     # Refitting of PLS model
#'     mypls <- PLS(
#'       xdata = x, ydata = y,
#'       selectedX = stab$selectedX,
#'       selectedY = stab$selectedY
#'     )
#'
#'     # Nonzero entries in weights are the same as in selectedX
#'     par(mfrow = c(2, 2))
#'     Heatmap(stab$selectedX,
#'       legend = FALSE
#'     )
#'     title("Selected in X")
#'     Heatmap(ifelse(mypls$Wmat != 0, yes = 1, no = 0),
#'       legend = FALSE
#'     )
#'     title("Nonzero entries in Wmat")
#'     Heatmap(stab$selectedY,
#'       legend = FALSE
#'     )
#'     title("Selected in Y")
#'     Heatmap(ifelse(mypls$Cmat != 0, yes = 1, no = 0),
#'       legend = FALSE
#'     )
#'     title("Nonzero entries in Cmat")
#'   }
#'
#'   # Multilevel PLS
#'   # Generating random design
#'   z <- rep(1:50, each = 4)
#'
#'   # Extracting the within-variability
#'   x_within <- mixOmics::withinVariation(X = x, design = cbind(z))
#'
#'   # Running PLS on within-variability
#'   mypls <- PLS(xdata = x_within, ydata = y, ncomp = 3)
#'
#'   par(oldpar)
#' }
#' }
#' @export
PLS <- function(xdata, ydata,
                selectedX = NULL, selectedY = NULL,
                family = "gaussian", ncomp = NULL, scale = TRUE) {
  # Checking mixOmics package is installed
  CheckPackageInstalled("mixOmics")

  # Checking arguments
  if (!family %in% c("gaussian")) {
    stop("Invalid input for argument 'family'. Only 'gaussian' family is supported.")
  }

  # Re-formatting ydata
  if (is.vector(ydata)) {
    ydata <- cbind(ydata)
  }
  if (is.null(colnames(ydata))) {
    colnames(ydata) <- paste0("outcome", 1:ncol(ydata))
  }

  # Checking consistency of arguments
  if (is.null(selectedX) & is.null(selectedY)) {
    if (is.null(ncomp)) {
      ncomp <- 1
    }
  } else {
    if (!is.null(selectedX)) {
      if (is.null(ncomp)) {
        ncomp <- ncol(selectedX)
      } else {
        if (!is.null(selectedY)) {
          if (ncol(selectedX) != ncol(selectedY)) {
            stop("Arguments 'selectedX' and 'selectedY' are not consistent. These matrices must have the same number of columns.")
          }
        }
        if (ncomp != ncol(selectedX)) {
          message(paste0(
            "Arguments 'ncomp' and 'selectedX' are not consistent. The algorithm is run with ncomp=",
            ncol(selectedX), "."
          ))
        }
      }
    } else {
      ncomp <- ncol(selectedY)
    }
  }

  # Defining the selection status for predictors
  if (is.null(selectedX)) {
    selectedX <- matrix(1, nrow = ncol(xdata), ncol = ncomp)
  }
  rownames(selectedX) <- colnames(xdata)
  colnames(selectedX) <- paste0("comp", 1:ncomp)

  # Defining the selection status for outcomes
  if (is.null(selectedY)) {
    selectedY <- matrix(1, nrow = ncol(ydata), ncol = ncomp)
  }
  rownames(selectedY) <- colnames(ydata)
  colnames(selectedY) <- paste0("comp", 1:ncomp)

  # Initialisation
  Emat <- scale(xdata, center = TRUE, scale = scale)
  Fmat <- scale(ydata, center = TRUE, scale = scale)
  meansX <- attr(Emat, "scaled:center")
  meansY <- attr(Fmat, "scaled:center")
  if (scale) {
    sigmaX <- attr(Emat, "scaled:scale")
    sigmaY <- attr(Fmat, "scaled:scale")
  } else {
    sigmaX <- rep(1, ncol(xdata))
    sigmaY <- rep(1, ncol(ydata))
  }

  # Initialisation of empty objects
  Wmat <- Pmat <- matrix(0, nrow = ncol(xdata), ncol = ncomp)
  rownames(Wmat) <- rownames(Pmat) <- colnames(xdata)
  Cmat <- matrix(0, nrow = ncol(ydata), ncol = ncomp)
  rownames(Cmat) <- colnames(ydata)
  Tmat <- matrix(NA, nrow = nrow(xdata), ncol = ncomp)
  Umat <- matrix(NA, nrow = nrow(ydata), ncol = ncomp)
  rownames(Tmat) <- rownames(Umat) <- rownames(xdata)
  colnames(Wmat) <- colnames(Pmat) <- colnames(Cmat) <- colnames(Tmat) <- colnames(Umat) <- paste0("comp", 1:ncomp)

  # Loop over the components
  for (comp in 1:ncomp) {
    # Extracting the stably selected predictors
    idsX <- rownames(selectedX)[which(selectedX[, comp] == 1)]
    if (length(idsX) > 0) {
      Emat_selected <- Emat[, idsX, drop = FALSE]
    } else {
      warning(paste0("No selected predictors for component ", comp, ". All predictors were included."))
      Emat_selected <- Emat
    }

    # Extracting the stably selected outcomes
    idsY <- rownames(selectedY)[which(selectedY[, comp] == 1)]
    if (length(idsY) > 0) {
      Fmat_selected <- Fmat[, idsY, drop = FALSE]
    } else {
      warning(paste0("No selected outcomes for component ", comp, ". All outcomes were included."))
      Fmat_selected <- Fmat
    }

    # Fitting PLS model (scaling done in initialisation and not required for further components)
    pls_h <- mixOmics::pls(
      X = Emat_selected,
      Y = Fmat_selected,
      mode = "regression",
      all.outputs = TRUE,
      multilevel = NULL,
      ncomp = 1,
      scale = FALSE
    )

    # Extracting the scores
    t_h <- pls_h$variates$X
    Tmat[, comp] <- t_h
    Umat[, comp] <- pls_h$variates$Y

    # Extracting the X-weights
    w_h <- pls_h$loadings$X
    Wmat[idsX, comp] <- w_h

    # Extracting the Y-weights
    Cmat[idsY, comp] <- pls_h$loadings$Y

    # Extracting the X-loadings
    Pmat[idsX, comp] <- pls_h$mat.c[, 1, drop = FALSE]

    # Deflation step
    Emat <- Emat - t_h %*% t(Pmat[, comp, drop = FALSE])
    d_h <- t(Fmat) %*% t_h * 1 / as.numeric(t(t_h) %*% t_h)
    Fmat <- Fmat - t_h %*% t(d_h)
  }

  # Transforming the weights
  Wstar <- base::zapsmall(Wmat %*% solve(t(Pmat) %*% Wmat))

  # Calculating matrices needed for predictions
  Qmat <- crossprod(scale(xdata, center = TRUE, scale = scale), Tmat)
  Rmat <- crossprod(scale(ydata, center = TRUE, scale = scale), Tmat)

  # Preparing outputs
  out <- list(
    Wmat = Wmat,
    Wstar = Wstar,
    Pmat = Pmat,
    Cmat = Cmat,
    Tmat = Tmat,
    Umat = Umat,
    Qmat = Qmat,
    Rmat = Rmat,
    meansX = meansX,
    sigmaX = sigmaX,
    meansY = meansY,
    sigmaY = sigmaY,
    methods = list(family = family, scale = scale),
    params = list(selectedX = selectedX, selectedY = selectedY)
  )

  return(out)
}


#' Partial Least Squares predictions
#'
#' Computes predicted values from a Partial Least Squares (PLS) model in
#' regression mode applied on \code{xdata}. This function is using the algorithm
#' implemented in \code{\link[mixOmics]{predict.pls}}.
#'
#' @inheritParams SelectionAlgo
#' @param model output of \code{\link{PLS}}.
#'
#' @return An array of predicted values.
#'
#' @seealso \code{\link{PLS}}
#'
#' @examples
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = c(5, 5, 5), family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'
#'   # PLS
#'   mypls <- PLS(xdata = x, ydata = y, ncomp = 3)
#'
#'   # Predicted values
#'   predicted <- PredictPLS(xdata = x, model = mypls)
#' }
#' @export
PredictPLS <- function(xdata, model) {
  # Extracting arguments
  family <- model$methods$family
  ncomp <- ncol(model$Wmat)

  # Checking arguments
  if (!family %in% c("gaussian")) {
    stop("Invalid input for argument 'family'. Only 'gaussian' family is supported.")
  }

  # Extracting relevant variables
  xdata <- xdata[, rownames(model$Wmat), drop = FALSE]

  # Re-scaling data
  if (ncol(xdata) == 1) {
    newdata <- matrix(sapply(xdata, FUN = function(x) {
      (x - model$meansX) / model$sigmaX
    }), ncol = 1)
  } else {
    newdata <- t(apply(xdata, 1, FUN = function(x) {
      (x - model$meansX) / model$sigmaX
    }))
  }

  # Initialising empty object
  out <- array(NA,
    dim = c(nrow(newdata), nrow(model$Rmat), ncomp),
    dimnames = list(rownames(newdata), rownames(model$Rmat), colnames(model$Rmat))
  )

  # Loop over the components
  for (comp in 1:ncomp) {
    # Computing predicted values
    Ypred <- newdata %*% model$Wmat[, 1:comp, drop = FALSE] %*% solve(t(model$Qmat[, 1:comp, drop = FALSE]) %*% model$Wmat[, 1:comp, drop = FALSE]) %*% t(model$Rmat)[1:comp, , drop = FALSE]

    # Re-scaling
    Ypred <- t(apply(Ypred, 1, FUN = function(y) {
      y * model$sigmaY + model$meansY
    }))

    # Storing in the output
    out[, , comp] <- Ypred
  }

  return(out)
}


#' Sparse Partial Least Squares
#'
#' Runs a sparse Partial Least Squares model using implementation from
#' \code{\link[sgPLS]{sgPLS-package}}. This function is not using stability.
#'
#' @inheritParams PLS
#' @param Lambda matrix of parameters controlling the number of selected
#'   predictors at current component, as defined by \code{ncomp}.
#' @param family type of PLS model. If \code{family="gaussian"}, a sparse PLS
#'   model as defined in \code{\link[sgPLS]{sPLS}} is run (for continuous
#'   outcomes). If \code{family="binomial"}, a PLS-DA model as defined in
#'   \code{\link[sgPLS]{sPLSda}} is run (for categorical outcomes).
#' @param keepX_previous number of selected predictors in previous components.
#'   Only used if \code{ncomp > 1}. The argument \code{keepX} in
#'   \code{\link[sgPLS]{sPLS}} is obtained by concatenating
#'   \code{keepX_previous} and \code{Lambda}.
#' @param keepY number of selected outcome variables. This argument is defined
#'   as in \code{\link[sgPLS]{sPLS}}. Only used if \code{family="gaussian"}.
#' @param scale logical indicating if the data should be scaled (i.e.
#'   transformed so that all variables have a standard deviation of one). Only
#'   used if \code{family="gaussian"}.
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
#' @family penalised dimensionality reduction functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @references \insertRef{sparsePLS}{sharp}
#'
#' @examples
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   ## Sparse PLS
#'
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 20, q = 3, family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'
#'   # Running sPLS with 2 X-variables and 1 Y-variable
#'   mypls <- SparsePLS(xdata = x, ydata = y, Lambda = 2, family = "gaussian", keepY = 1)
#'
#'
#'   ## Sparse PLS-DA
#'
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 20, family = "binomial")
#'
#'   # Running sPLS-DA with 2 X-variables and 1 Y-variable
#'   mypls <- SparsePLS(xdata = simul$xdata, ydata = simul$ydata, Lambda = 2, family = "binomial")
#' }
#' @export
SparsePLS <- function(xdata, ydata, Lambda, family = "gaussian",
                      ncomp = 1, scale = TRUE,
                      keepX_previous = NULL, keepY = NULL, ...) {
  # Checking sgPLS package is installed
  CheckPackageInstalled("sgPLS")

  # Storing extra arguments
  extra_args <- list(...)

  # Checking arguments
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
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
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::sPLS)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "scale", "keepX", "keepY")]

      # Running PLS model
      mymodel <- do.call(sgPLS::sPLS, args = c(
        list(
          X = xdata, Y = ydata,
          ncomp = ncomp, scale = scale,
          keepX = nvarx, keepY = keepY
        ),
        tmp_extra_args
      ))
    } else {
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::sPLSda)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "keepX")]

      # Running PLS model
      mymodel <- do.call(sgPLS::sPLSda, args = c(
        list(
          X = xdata, Y = as.vector(ydata),
          ncomp = ncomp,
          keepX = nvarx
        ),
        tmp_extra_args
      ))
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
#' @family penalised dimensionality reduction functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @references \insertRef{sparsegroupPLS}{sharp}
#'
#' @examples
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   ## Sparse group PLS
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 30, q = 3, family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'
#'   # Running sgPLS with 1 group and sparsity of 0.5
#'   mypls <- SparseGroupPLS(
#'     xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'     group_x = c(10, 15, 5), alpha.x = 0.5
#'   )
#'
#'   # Running sgPLS with groups on outcomes
#'   mypls <- SparseGroupPLS(
#'     xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'     group_x = c(10, 15, 5), alpha.x = 0.5,
#'     group_y = c(2, 1), keepY = 1, alpha.y = 0.9
#'   )
#'
#'   ## Sparse group PLS-DA
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 50, family = "binomial")
#'
#'   # Running sgPLS-DA with 1 group and sparsity of 0.9
#'   mypls <- SparseGroupPLS(
#'     xdata = simul$xdata, ydata = simul$ydata, Lambda = 1, family = "binomial",
#'     group_x = c(10, 15, 25), alpha.x = 0.9
#'   )
#' }
#' @export
SparseGroupPLS <- function(xdata, ydata, family = "gaussian", group_x, group_y = NULL,
                           Lambda, alpha.x, alpha.y = NULL,
                           keepX_previous = NULL, keepY = NULL,
                           ncomp = 1, scale = TRUE, ...) {
  # Checking sgPLS package is installed
  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Checking arguments
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  # Checking groups
  if (sum(group_x) != ncol(xdata)) {
    stop("Invalid input for argument 'group_x'. The sum of the number of variables per group must be equal to the total number of variables in 'xdata'.")
  }
  if (!is.null(group_y)) {
    if (sum(group_y) != ncol(ydata)) {
      stop("Invalid input for argument 'group_y'. The sum of the number of variables per group must be equal to the total number of variables in 'ydata'.")
    }
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
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::sgPLSda)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "ind.block.x", "keepX", "alpha.x")]

      # Running PLS model
      mymodel <- do.call(sgPLS::sgPLSda, args = c(
        list(
          X = xdata, Y = as.vector(ydata),
          ncomp = ncomp,
          ind.block.x = ind.block.x,
          keepX = nvarx, alpha.x = alpha.x
        ),
        tmp_extra_args
      )) # no sparsity in Y
    } else {
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::sgPLS)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "scale", "ind.block.x", "ind.block.y", "keepX", "keepY", "alpha.x", "alpha.y")]

      # Running PLS model
      mymodel <- do.call(sgPLS::sgPLS, args = c(
        list(
          X = xdata, Y = ydata,
          ncomp = ncomp, scale = scale,
          ind.block.x = ind.block.x, ind.block.y = ind.block.y,
          keepX = nvarx, keepY = keepY,
          alpha.x = alpha.x, alpha.y = alpha.y
        ),
        tmp_extra_args
      ))
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
#' @family penalised dimensionality reduction functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @references \insertRef{sparsegroupPLS}{sharp}
#'
#' @examples
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   ## Group PLS
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 50, q = 3, family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'
#'   # Running gPLS with 1 group and sparsity of 0.5
#'   mypls <- GroupPLS(
#'     xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'     group_x = c(10, 15, 25),
#'   )
#'
#'   # Running gPLS with groups on outcomes
#'   mypls <- GroupPLS(
#'     xdata = x, ydata = y, Lambda = 1, family = "gaussian",
#'     group_x = c(10, 15, 25),
#'     group_y = c(2, 1), keepY = 1
#'   )
#' }
#' @export
GroupPLS <- function(xdata, ydata, family = "gaussian", group_x, group_y = NULL,
                     Lambda, keepX_previous = NULL, keepY = NULL,
                     ncomp = 1, scale = TRUE, ...) {
  # Checking sgPLS package is installed
  if (!requireNamespace("sgPLS")) {
    stop("This function requires the 'sgPLS' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Checking arguments
  if (!family %in% c("binomial", "gaussian")) {
    stop("Invalid input for argument 'family'. For PLS models, argument 'family' must be 'gaussian' or 'binomial'.")
  }

  # Checking groups
  if (sum(group_x) != ncol(xdata)) {
    stop("Invalid input for argument 'group_x'. The sum of the number of variables per group must be equal to the total number of variables in 'xdata'.")
  }
  if (!is.null(group_y)) {
    if (sum(group_y) != ncol(ydata)) {
      stop("Invalid input for argument 'group_y'. The sum of the number of variables per group must be equal to the total number of variables in 'ydata'.")
    }
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
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::gPLSda)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "ind.block.x", "keepX")]

      # Running PLS model
      mymodel <- do.call(sgPLS::gPLSda, args = c(
        list(
          X = xdata, Y = as.vector(ydata),
          ncomp = ncomp,
          ind.block.x = ind.block.x,
          keepX = nvarx
        ),
        tmp_extra_args
      )) # no sparsity in Y
    } else {
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sgPLS::gPLSda)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "Y", "ncomp", "ind.block.x", "keepX")]

      # Running PLS model
      mymodel <- do.call(sgPLS::gPLS, args = c(
        list(
          X = xdata, Y = ydata,
          ncomp = ncomp, scale = scale,
          ind.block.x = ind.block.x, ind.block.y = ind.block.y,
          keepX = nvarx, keepY = keepY
        ),
        tmp_extra_args
      ))
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


#' Sparse Principal Component Analysis
#'
#' Runs a sparse Principal Component Analysis model using implementation from
#' \code{\link[elasticnet]{spca}} (if \code{algo="sPCA"}) or
#' \code{\link[mixOmics]{spca}} (if \code{algo="rSVD"}). This function is not
#' using stability.
#'
#' @inheritParams PLS
#' @param xdata data matrix with observations as rows and variables as columns.
#' @param Lambda matrix of parameters controlling the number of selected
#'   variables at current component, as defined by \code{ncomp}.
#' @param keepX_previous number of selected predictors in previous components.
#'   Only used if \code{ncomp > 1}.
#' @param algorithm character string indicating the name of the algorithm to use for
#'   sparse PCA. Possible values are: "sPCA" (for the algorithm proposed by Zou,
#'   Hastie and Tibshirani and implemented in \code{\link[elasticnet]{spca}}) or
#'   "rSVD" (for the regularised SVD approach proposed by Shen and Huang and
#'   implemented in \code{\link[mixOmics]{spca}}).
#' @param ... additional arguments to be passed to
#'   \code{\link[elasticnet]{spca}} (if \code{algorithm="sPCA"}) or
#'   \code{\link[mixOmics]{spca}} (if \code{algorithm="rSVD"}).
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors (starting
#'   with "X") or outcomes (starting with "Y") variables for different
#'   components (denoted by "PC").}
#'
#' @family penalised dimensionality reduction functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @references \insertRef{sparsePCA}{sharp}
#'
#'   \insertRef{sparsePCASVD}{sharp}
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' x <- simul$xdata
#'
#' # Sparse PCA (by Zou, Hastie, Tibshirani)
#' if (requireNamespace("elasticnet", quietly = TRUE)) {
#'   mypca <- SparsePCA(
#'     xdata = x, ncomp = 2,
#'     Lambda = c(1, 2), keepX_previous = 10, algorithm = "sPCA"
#'   )
#' }
#'
#' # Sparse PCA (by Shen and Huang)
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   mypca <- SparsePCA(
#'     xdata = x, ncomp = 2,
#'     Lambda = c(1, 2), keepX_previous = 10, algorithm = "rSVD"
#'   )
#' }
#' @export
SparsePCA <- function(xdata, Lambda,
                      ncomp = 1, scale = TRUE,
                      keepX_previous = NULL,
                      algorithm = "sPCA", ...) {
  # Checking input value for algorithm
  if (!algorithm %in% c("sPCA", "rSVD")) {
    stop("Invalid input for argument 'algorithm'. Possible values are: 'sPCA' or 'rSVD'.")
  }

  # Checking mixOmics package is installed
  if (algorithm == "rSVD") {
    if (!requireNamespace("mixOmics")) {
      stop("This function requires the 'mixOmics' package.")
    }
  }

  # Checking elasticnet package is installed
  if (algorithm == "sPCA") {
    if (!requireNamespace("elasticnet")) {
      stop("This function requires the 'elasticnet' package.")
    }
  }

  # Storing extra arguments
  extra_args <- list(...)

  # All X variables are kept in previous components by default
  if (is.null(keepX_previous)) {
    keepX_previous <- rep(ncol(xdata), ncomp - 1)
  }

  # Initialising the current set of loadings coefficients
  beta <- matrix(NA, nrow = length(Lambda), ncol = ncol(xdata))
  rownames(beta) <- paste0("s", 0:(length(Lambda) - 1))
  colnames(beta) <- colnames(xdata)

  # # Initialising the full set of loadings coefficients
  beta_full <- matrix(NA, nrow = length(Lambda), ncol = ncol(xdata) * ncomp)
  rownames(beta_full) <- paste0("s", 0:(length(Lambda) - 1))
  colnames(beta_full) <- c(paste0(paste0("X_", colnames(xdata), "_PC"), rep(1:ncomp, each = ncol(xdata))))

  # Loop over all parameters (number of selected variables in X in current component)
  for (k in 1:length(Lambda)) {
    # Number of selected variables per component in X
    nvarx <- c(keepX_previous, Lambda[k])

    if (algorithm == "rSVD") {
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = mixOmics::spca)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "ncomp", "keepX", "scale")]

      # Running PLS model
      mymodel <- do.call(mixOmics::spca, args = c(
        list(
          X = xdata,
          ncomp = ncomp,
          scale = scale,
          keepX = nvarx
        ),
        tmp_extra_args
      ))

      # Extracting X and Y loadings
      Xloadings <- mymodel$loadings$X
    }

    if (algorithm == "sPCA") {
      # Extracting relevant extra arguments
      tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = elasticnet::spca)
      tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "type", "K", "para", "sparse")]

      # Running PLS model
      if (scale) {
        mymodel <- do.call(elasticnet::spca, args = c(list(
          x = stats::cor(xdata), type = "Gram",
          K = ncomp, para = nvarx, sparse = "varnum"
        ), tmp_extra_args))
      } else {
        mymodel <- do.call(elasticnet::spca, args = c(list(
          x = stats::cov(xdata), type = "Gram",
          K = ncomp, para = nvarx, sparse = "varnum"
        ), tmp_extra_args))
      }

      # Extracting X and Y loadings
      Xloadings <- mymodel$loadings
    }

    # Exporting the set of loadings to look at for selection with current component
    beta[k, ] <- Xloadings[, ncomp, drop = FALSE]

    # Exporting all loadings coefficients
    beta_full[k, ] <- as.vector(Xloadings)
  }

  beta <- ifelse(beta != 0, yes = 1, no = 0)

  return(list(selected = beta, beta_full = beta_full))
}

#' Calibrated parameter indices
#'
#' Extracts the indices of calibrated parameters with respect to the grids provided
#' in \code{Lambda} and \code{pi_list} in \code{stability}.
#'
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}. If \code{stability=NULL}, \code{S} must be
#'   provided.
#' @param S matrix of stability scores obtained with different combinations of
#'   parameters where rows correspond to different values of the parameter
#'   controlling the level of sparsity in the underlying feature selection
#'   algorithm and columns correspond to different values of the threshold in
#'   selection proportions. If \code{S=NULL}, argument \code{stability} must be
#'   provided.
#'
#' @return A matrix of parameter indices. In multi-block graphical modelling, rows
#'   correspond to different blocks.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Extracting IDs of calibrated parameters
#' ids <- ArgmaxId(stab)
#' stab$Lambda[ids[1], 1]
#' stab$params$pi_list[ids[2]]
#'
#' # Link with Argmax() function
#' args <- Argmax(stab)
#' }
#'
#' @export
ArgmaxId <- function(stability = NULL, S = NULL) {
  if ((is.null(stability)) & (is.null(S))) {
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }
  if (is.null(S)) {
    argmax_id <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
    if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
      id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
      argmax_id[, 1] <- rep(id, nrow(argmax_id))
      for (block_id in 1:ncol(stability$Lambda)) {
        if (!is.na(stability$P[id, block_id])) {
          argmax_id[block_id, 2] <- which(stability$params$pi_list == stability$P[id, block_id])
        }
      }
    } else {
      for (block_id in 1:ncol(stability$Lambda)) {
        if (ncol(stability$Lambda) == 1) {
          myS <- stability$S
        } else {
          myS <- stability$S[, block_id, drop = FALSE]
        }
        myS[is.na(myS)] <- 0
        myid <- which.max(myS[, 1])
        argmax_id[block_id, ] <- c(myid, which(stability$params$pi_list == stability$P[myid, block_id]))
      }
    }
  } else {
    argmax_id <- matrix(NA, nrow = 1, ncol = 2)
    myS <- apply(S, 1, max, na.rm = TRUE)
    myS[is.na(myS)] <- 0
    myid <- which.max(myS)
    argmax_id[1, ] <- c(myid, max(which(S[myid, ] == myS[myid])))
  }
  colnames(argmax_id) <- c("lambda_id", "pi_id")
  return(argmax_id)
}


#' Calibrated parameters
#'
#' Extracts calibrated parameter values in stability selection.
#'
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}.
#'
#' @return A matrix of parameter values. In multi-block graphical modelling,
#'   rows correspond to different blocks.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Extracting calibrated parameters
#' args <- Argmax(stab)
#' }
#'
#' @export
Argmax <- function(stability) {
  argmax <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
  if (is.null(stability$params$lambda_other_blocks) & (length(stability$params$pk) > 1)) {
    id <- which.max(apply(stability$S, 1, sum, na.rm = TRUE))
    argmax[, 1] <- stability$Lambda[id, ]
    argmax[, 2] <- stability$P[id, ]
  } else {
    for (block_id in 1:ncol(stability$Lambda)) {
      if (ncol(stability$Lambda) == 1) {
        myS <- stability$S
      } else {
        myS <- stability$S[, block_id, drop = FALSE]
      }
      myS[is.na(myS)] <- 0
      myid <- which.max(myS[, 1])
      argmax[block_id, ] <- c(stability$Lambda[myid, block_id], stability$P[myid, block_id])
    }
  }
  colnames(argmax) <- c("lambda", "pi")
  return(argmax)
}


#' Calibrated adjacency matrix
#'
#' Builds the adjacency matrix of the (calibrated) stability selection graphical
#' model.
#'
#' @param stability output of \code{\link{GraphicalModel}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated adjacency matrix is returned.
#'
#' @return A binary and symmetric adjacency matrix encoding an undirected graph
#'   with no self-loops.
#'
#' @family calibration functions
#' @seealso \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Calibrated adjacency matrix
#' A <- Adjacency(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(20, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' A <- Adjacency(stab, argmax_id = myids)
#' }
#'
#' @export
Adjacency <- function(stability, argmax_id = NULL) {
  A <- matrix(0, ncol = ncol(stability$selprop), nrow = nrow(stability$selprop))
  bigblocks <- BlockMatrix(stability$params$pk)
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
    argmax <- Argmax(stability)
  } else {
    argmax <- NULL
    for (block_id in 1:ncol(stability$Lambda)) {
      argmax <- rbind(argmax, c(
        stability$Lambda[argmax_id[block_id, 1], block_id],
        stability$params$pi_list[argmax_id[block_id, 2]]
      ))
    }
  }
  for (block_id in 1:ncol(stability$Lambda)) {
    A_block <- ifelse(stability$selprop[, , argmax_id[block_id, 1]] >= argmax[block_id, 2], 1, 0)
    A_block[lower.tri(A_block)] <- 0
    A_block <- A_block + t(A_block) # for symmetry
    if (length(stability$params$pk) > 1) {
      A_block[bigblocks != block_id] <- 0
    }
    A <- A + A_block
  }
  A[is.na(A)] <- 0
  return(A)
}


#' Set of stably selected variables
#'
#' Builds the (calibrated) set of stably selected variables.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}}.
#'
#' @return A binary vector encoding the selection status of the variables
#'   (\code{1} if selected, \code{0} otherwise).
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$X, ydata = simul$Y)
#'
#' # Calibrated set
#' selected <- SelectedVariables(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(80, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' selected <- SelectedVariables(stab, argmax_id = myids)
#' }
#' @export
SelectedVariables <- function(stability, argmax_id = NULL) {
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
  }
  stability_selected <- ifelse(stability$selprop[argmax_id[1], ] >= stability$params$pi_list[argmax_id[2]], 1, 0)
  return(stability_selected)
}


#' Stable cluster membership
#'
#' Builds the (calibrated) stable clusters as connected components of the graph
#' defined from stable co-membership.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{Clustering}}.
#'
#' @return A vector encoding the cluster membership.
#'
#' @family calibration functions
#' @seealso \code{\link{Clustering}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 50, n = 10)
#'
#' # Stability selection
#' stab <- Clustering(xdata = simul$data, Lambda = 1:nrow(simul$data))
#'
#' # Stable cluster membership
#' mymembership <- Clusters(stab)
#' }
#' @export
Clusters <- function(stability, argmax_id = NULL) {
  # Computing stable co-membership matrix
  adjacency <- Adjacency(stability = stability, argmax_id = argmax_id)

  # Extracting stable connected components
  mymembership <- igraph::components(Graph(adjacency, satellites = TRUE))$membership

  return(mymembership)
}


#' Selection proportions
#'
#' Extracts the selection proportions of the (calibrated) stability selection
#' model.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}. If \code{stability=NULL}, \code{S} must be
#'   provided.
#'
#' @return A symmetric matrix (graphical model) or vector (variable selection)
#'   of selection proportions.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' ## Variable selection
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$X, ydata = simul$Y)
#'
#' # Calibrated selection proportions
#' prop <- SelectionProportions(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(80, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' prop <- SelectionProportions(stab, argmax_id = myids)
#'
#' ## Graphical model
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Calibrated adjacency matrix
#' prop <- SelectionProportions(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(20, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' prop <- SelectionProportions(stab, argmax_id = myids)
#' }
#'
#' @export
SelectionProportions <- function(stability, argmax_id = NULL) {
  if ("Sequential_template" %in% names(stability$params)) {
    out <- SelectionProportionsGraphical(stability = stability, argmax_id = argmax_id)
  } else {
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  return(out)
}


#' Selection proportions (graphical model)
#'
#' Extracts the selection proportions of the (calibrated) stability selection
#' model.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{GraphicalModel}}.
#'
#' @return A symmetric matrix.
#'
#' @keywords internal
SelectionProportionsGraphical <- function(stability, argmax_id = NULL) {
  A <- matrix(0, ncol = ncol(stability$selprop), nrow = nrow(stability$selprop))
  bigblocks <- BlockMatrix(stability$params$pk)
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
    argmax <- Argmax(stability)
  } else {
    argmax <- NULL
    for (block_id in 1:ncol(stability$Lambda)) {
      argmax <- rbind(argmax, c(
        stability$Lambda[argmax_id[block_id, 1], ],
        stability$params$pi_list[argmax_id[block_id, 2]]
      ))
    }
  }
  for (block_id in 1:ncol(stability$Lambda)) {
    A_block <- stability$selprop[, , argmax_id[block_id, 1]]
    A_block[lower.tri(A_block)] <- 0
    A_block <- A_block + t(A_block) # for symmetry
    if (length(stability$params$pk) > 1) {
      A_block[bigblocks != block_id] <- 0
    }
    A <- A + A_block
  }
  return(A)
}


#' Selection proportions (variable selection)
#'
#' Extracts the selection proportions of the (calibrated) stability selection
#' model.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}}.
#'
#' @return A vector of selection proportions.
#'
#' @keywords internal
SelectionProportionsRegression <- function(stability, argmax_id = NULL) {
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
    argmax <- Argmax(stability)
  }
  m <- stability$selprop[argmax_id[1], ]
  calibrated_pi <- stability$params$pi_list[argmax_id[2]]
  return(m)
}


#' Model coefficients
#'
#' Extracts the coefficients of visited models at different resampling
#' iterations of a stability selection run. This function can be applied to the
#' output of \code{\link{VariableSelection}}.
#'
#' @param stability output of \code{\link{VariableSelection}}.
#' @param side character string indicating if coefficients of the predictor
#'   (\code{side="X"}) or outcome (\code{side="Y"}) coefficients should be
#'   returned. Option \code{side="Y"} is only applicable to PLS models.
#' @param comp component ID. Only applicable to PLS models.
#' @param iterations vector of iteration IDs. If \code{iterations=NULL}, the
#'   coefficients of all visited models are returned. This number must be
#'   smaller than the number of iterations \code{K} used for the stability
#'   selection run.
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$X, ydata = simul$Y, family = "gaussian")
#'
#' # Coefficients of visited models
#' coefs <- Coefficients(stab)
#' dim(coefs)
#'
#' # Coefficients of the first fitted model
#' coefs <- Coefficients(stab, iterations = 1)
#' dim(coefs)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = simul$X, ydata = simul$Y,
#'   implementation = SparsePLS, family = "gaussian"
#' )
#'
#' # Coefficients of visited models
#' coefs <- Coefficients(stab, side = "Y", )
#' dim(coefs)
#' }
#' @export
Coefficients <- function(stability, side = "X", comp = 1, iterations = NULL) {
  if (is.null(iterations)) {
    iterations <- seq(1, dim(stability$Beta)[3])
  } else {
    iterations <- iterations[iterations <= dim(stability$Beta)[3]]
  }
  if (length(iterations) == 0) {
    stop("Invalid input for argument 'iterations'. This argument must be a number smaller than the number of iterations used for the stability selection run.")
  }
  if (ncol(stability$Beta) == stability$params$pk) {
    return(stability$Beta[, , iterations, drop = FALSE])
  } else {
    if (!side %in% c("X", "Y")) {
      warning("Invalid input for argument 'side'. The default value ('X') was used.")
      side <- "X"
    }
    side_id <- grepl(paste0(side, "_"), colnames(stability$Beta))
    comp_id <- grepl(paste0("_PC", comp), colnames(stability$Beta))
    return(stability$Beta[, side_id & comp_id, iterations, drop = FALSE])
  }
}

#' Stable results
#'
#' Extracts stable results for stability selection or consensus clustering.
#'
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{BiSelection}}, \code{\link{GraphicalModel}} or
#'   \code{\link{Clustering}}.
#' @param argmax_id optional indices of hyper-parameters. If
#'   \code{argmax_id=NULL}, the calibrated hyper-parameters are used.
#' @param linkage character string indicating the type of linkage used in
#'   hierarchical clustering to define the stable clusters. Possible values
#'   include \code{"complete"}, \code{"single"} and \code{"average"} (see
#'   argument \code{"method"} in \code{\link[stats]{hclust}} for a full list).
#'
#' @return A binary vector or matrix encoding the selection status (\code{1} if
#'   selected, \code{0} otherwise) in stability selection or stable cluster
#'   membership in consensus clustering.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}},
#'   \code{\link{GraphicalModel}}, \code{\link{Clustering}}
#'
#' @examples
#' \donttest{
#' # Variable selection
#' set.seed(1)
#' simul <- SimulateRegression(pk = 20)
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' SelectedVariables(stab)
#' Stable(stab)
#'
#' # Graphical model
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 10)
#' stab <- GraphicalModel(xdata = simul$data)
#' Adjacency(stab)
#' Stable(stab)
#'
#' # Clustering
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(30, 30, 30),
#'   nu_xc = 1
#' )
#' stab <- Clustering(xdata = simul$data)
#' Clusters(stab)
#' Stable(stab)
#' }
#' @export
Stable <- function(stability, argmax_id = NULL, linkage = "complete") {
  if (inherits(stability, c("variable_selection", "bi_selection"))) {
    out <- SelectedVariables(stability = stability, argmax_id = argmax_id)
  }

  if (inherits(stability, c("graphical_model"))) {
    out <- Adjacency(stability = stability, argmax_id = argmax_id)
  }

  if (inherits(stability, c("clustering"))) {
    out <- Clusters(stability = stability, linkage = linkage, argmax_id = argmax_id)
  }

  return(out)
}


#' @rdname Stable
#' @export
SelectedVariables <- function(stability, argmax_id = NULL) {
  if (!inherits(stability, c("variable_selection", "bi_selection"))) {
    stop("Invalid input for argument 'stability'. This function only applies to outputs from VariableSelection() or BiSelection().")
  }

  if (inherits(stability, "variable_selection")) {
    if (is.null(argmax_id)) {
      argmax_id <- ArgmaxId(stability)
      argmax <- Argmax(stability)
    } else {
      argmax <- c(NA, stability$params$pi_list[argmax_id[2]])
    }
    stability_selected <- ifelse(stability$selprop[argmax_id[1], ] >= argmax[2],
      yes = 1, no = 0
    )
  }

  if (inherits(stability, "bi_selection")) {
    if (is.null(argmax_id)) {
      stability_selected <- stability$selectedX
    } else {
      stop("Invalid input for argument 'argmax_id'. Arbitrary choice of parameters is not supported for the output of BiSelection().")
    }
  }

  return(stability_selected)
}


#' @rdname Stable
#' @export
Adjacency <- function(stability, argmax_id = NULL) {
  if (!inherits(stability, c("graphical_model", "bi_selection", "clustering"))) {
    stop("Invalid input for argument 'stability'. This function only applies to outputs from GraphicalModel(), BiSelection() or Clustering().")
  }

  if (inherits(stability, "bi_selection")) {
    if ("selectedY" %in% names(stability)) {
      A <- Square(t(rbind(stability$selectedX, stability$selectedY)))
    } else {
      A <- Square(stability$selectedX)
    }
  }

  if (inherits(stability, "graphical_model")) {
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
      if (length(stability$params$pk) > 1) {
        A_block[bigblocks != block_id] <- 0
      }
      A <- A + A_block
    }
  }

  if (inherits(stability, "clustering")) {
    A <- CoMembership(Clusters(stability = stability))
  }

  A[is.na(A)] <- 0
  return(A)
}


#' @rdname Stable
#' @export
Clusters <- function(stability, linkage = "complete", argmax_id = NULL) {
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability = stability)
  }

  # Calibrated consensus matrix
  coprop <- ConsensusMatrix(stability = stability, argmax_id = argmax_id[1])

  # Extracting stable clusters from hierarchical clustering
  shclust <- stats::hclust(stats::as.dist(1 - coprop), method = linkage)
  mymembership <- stats::cutree(shclust, k = ceiling(stability$nc[argmax_id[1], 1]))

  # Splitting noise clusters
  if (any(!is.na(stability$bignoise))) {
    noise_prop <- stability$bignoise[, argmax_id[1]]
    if (any(noise_prop >= 0.5)) {
      ids <- which(noise_prop >= 0.5)
      mymembership[ids] <- max(mymembership) + seq(1, length(ids))
      mymembership <- mymembership - min(mymembership) + 1
    }
  }

  return(mymembership)
}


#' Selection/co-membership proportions
#'
#' Extracts selection proportions (for stability selection) or co-membership
#' proportions (for consensus clustering).
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}}, \code{\link{BiSelection}}, or
#'   \code{\link{Clustering}}.
#'
#' @return A vector or matrix of proportions.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}, \code{\link{Clustering}}
#'
#' @examples
#' \donttest{
#' # Stability selection
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' SelectionProportions(stab)
#'
#' # Consensus clustering
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(30, 30, 30), nu_xc = 1, ev_xc = 0.5
#' )
#' stab <- Clustering(xdata = simul$data)
#' ConsensusMatrix(stab)
#' }
#'
#' @export
SelectionProportions <- function(stability, argmax_id = NULL) {
  out <- NULL

  if (inherits(stability, "graphical_model")) {
    out <- SelectionProportionsGraphical(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "variable_selection")) {
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "clustering")) {
    argmax_id <- ArgmaxId(stability)
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (inherits(stability, "bi_selection")) {
    out <- stability$selpropX
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
  }
  m <- stability$selprop[argmax_id[1], ]
  return(m)
}


#' @rdname SelectionProportions
#' @export
ConsensusMatrix <- function(stability, argmax_id = NULL) {
  if (inherits(stability, "clustering")) {
    if (is.null(argmax_id)) {
      argmax_id <- ArgmaxId(stability = stability)
    }
    mat <- stability$coprop[, , argmax_id[1]]
  } else {
    stop("Invalid input for argument 'stability'. Only applicable to an object of class 'clustering', i.e. to the output of Clustering().")
  }

  return(mat)
}

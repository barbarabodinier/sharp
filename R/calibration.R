#' Calibrated parameter indices
#'
#' Extracts the indices of calibrated parameters with respect to the grids
#' provided in \code{Lambda} and \code{pi_list} in \code{stability}.
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
#' @return A matrix of parameter indices. For multi-block graphical models, rows
#'   correspond to different blocks.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
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
#' # Alternative formulation
#' ids2 <- ArgmaxId(S = stab$S_2d)
#'
#' # Link with Argmax() function
#' args <- Argmax(stab)
#' @export
ArgmaxId <- function(stability = NULL, S = NULL) {
  if ((is.null(stability)) & (is.null(S))) {
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }

  clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)

  if (is.null(S)) {
    if (clustering) {
      # If multiple solutions, prioritising many clusters over small number of features
      Sc <- round(stability$Sc, digits = 4)
      argmax_id <- which(Sc == max(Sc, na.rm = TRUE))
      argmax_id <- argmax_id[which(stability$nc[argmax_id] == max(stability$nc[argmax_id]))]
      # argmax_id <- argmax_id[which(stability$Q[argmax_id] == min(stability$Q[argmax_id]))]
      argmax_id <- cbind(min(argmax_id))

      # if (any(!is.na(stability$S_2d[argmax_id, ]))) {
      #   argmax_id <- matrix(c(argmax_id, which.max(stability$S_2d[argmax_id, ])), ncol = 2)
      # } else {
      #   argmax_id <- matrix(c(argmax_id, NA), ncol = 2)
      # }
    } else {
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
    }
  } else {
    argmax_id <- matrix(NA, nrow = 1, ncol = 2)
    myS <- apply(S, 1, max, na.rm = TRUE)
    myS[is.na(myS)] <- 0
    myid <- which.max(myS)
    argmax_id[1, ] <- c(myid, max(which(S[myid, ] == myS[myid])))
  }
  if (clustering) {
    colnames(argmax_id) <- c("row_id")
  } else {
    colnames(argmax_id) <- c("lambda_id", "pi_id")
  }
  return(argmax_id)
}


#' Calibrated parameters
#'
#' Extracts calibrated parameter values in stability selection.
#'
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{BiSelection}} or \code{\link{GraphicalModel}}.
#'
#' @return A matrix of parameter values. If applied to the output of
#'   \code{\link{VariableSelection}} or \code{\link{GraphicalModel}}, the first
#'   column (\code{lambda}) denotes the calibrated hyper-parameter of the
#'   underlying algorithm. The second column (\code{pi}) is the calibrated
#'   threshold in selection/co-membership proportions. For multi-block graphical
#'   models, rows correspond to different blocks. If applied to the output of
#'   \code{\link{BiSelection}}, all columns are named as in object
#'   \code{summary}.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}
#'
#' @examples
#' ## Graphical modelling
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Extracting calibrated parameters
#' Argmax(stab)
#' @export
Argmax <- function(stability) {
  clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)
  if (inherits(stability, "bi_selection")) {
    argmax <- stability$summary
    argmax <- argmax[, colnames(argmax) != "S", drop = FALSE]
  } else {
    if (clustering) {
      id <- ArgmaxId(stability = stability)
      argmax <- matrix(c(
        stability$nc[id[1], 1],
        stability$Lambda[id[1], 1]
      ),
      ncol = 2
      )
    } else {
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
    }
    if (clustering) {
      colnames(argmax) <- c("nc", "lambda")
    } else {
      colnames(argmax) <- c("lambda", "pi")
    }
  }

  return(argmax)
}


#' Calibrated adjacency matrix
#'
#' Extracts the adjacency matrix of the (calibrated) stability selection
#' graphical model.
#'
#' @param stability output of \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}, or \code{\link{Clustering}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated model is used.
#'
#' @return A binary and symmetric adjacency matrix encoding an undirected graph
#'   with no self-loops.
#'
#' @family calibration functions
#' @seealso \code{\link{GraphicalModel}}
#'
#' @examples
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


#' Stably selected variables
#'
#' Extracts the (calibrated) set of stably selected variables.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}}, or
#'   \code{\link{BiSelection}}.
#'
#' @return A binary vector encoding the selection status of the variables
#'   (\code{1} if selected, \code{0} otherwise).
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#'
#' # Calibrated set
#' selected <- SelectedVariables(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(50, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' selected <- SelectedVariables(stab, argmax_id = myids)
#' }
#' @export
SelectedVariables <- function(stability, argmax_id = NULL) {
  if (!inherits(stability, c("clustering", "variable_selection", "bi_selection"))) {
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


#' Selection proportions
#'
#' Extracts the selection (or co-membership) proportions of the (calibrated)
#' model.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}}, or \code{\link{BiSelection}}.
#'
#' @return A symmetric matrix (graphical model) or vector (variable selection)
#'   of selection proportions.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}
#'
#' @examples
#' \donttest{
#' ## Variable selection
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
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
#'
#' ## Graphical model
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Calibrated matrix of selection proportions
#' prop <- SelectionProportions(stab)
#'
#' # User-defined parameters
#' myids <- matrix(c(20, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' prop <- SelectionProportions(stab, argmax_id = myids)
#'
#'
#' ## Dimensionality reduction
#'
#' # Data simulation (continuous outcomes)
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#' x <- simul$xdata
#' y <- simul$ydata
#'
#' # Sparse PLS
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#'
#' # Calibrated selection proportions per component
#' prop <- SelectionProportions(stab)
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
#' @keywords internal
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


#' Summarised coefficients conditionally on selection
#'
#' Computes descriptive statistics (defined by \code{FUN}) for coefficients of
#' the (calibrated) models conditionally on selection across resampling
#' iterations.
#'
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{BiSelection}}.
#' @param side character string indicating if coefficients of predictors
#'   (\code{side="X"}) or outcomes (\code{side="Y"}) should be returned. Only
#'   applicable to PLS models.
#' @param comp component ID. Only applicable to PLS models.
#' @param lambda_id parameter ID with respect to the grid \code{Lambda}. If
#'   \code{NULL}, aggregated coefficients across the models run with the
#'   calibrated parameter are returned.
#' @param FUN function to use to aggregate coefficients of visited models over
#'   resampling iterations. Recommended functions include
#'   \code{\link[stats]{median}} or \code{\link[base]{mean}}.
#' @param ... additional arguments to be passed to \code{FUN}.
#'
#' @return A matrix of summarised coefficients conditionally on selection across
#'   resampling iterations. Missing values (\code{NA}) are returned for
#'   variables that are never selected.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}},
#'   \code{\link{Refit}}
#'
#' @examples
#' \donttest{
#' # Example with univariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' median_betas <- AggregatedEffects(stab)
#'
#' # Comparison with refitted model
#' refitted <- Refit(xdata = simul$xdata, ydata = simul$ydata, stability = stab)
#' refitted_betas <- refitted$coefficients[-1]
#' plot(median_betas[names(refitted_betas), ], refitted_betas,
#'   panel.first = abline(0, 1, lty = 2)
#' )
#'
#' # Extracting mean betas conditionally on selection
#' mean_betas <- AggregatedEffects(stab, FUN = mean)
#' plot(median_betas, mean_betas)
#'
#' # Regression with multivariate outcomes
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, q = 2, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "mgaussian")
#' median_betas <- AggregatedEffects(stab)
#' dim(median_betas)
#'
#' # Sparse PLS with multivariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = 15, q = 3, family = "gaussian")
#' x <- simul$xdata
#' y <- simul$ydata
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#' median_betas <- AggregatedEffects(stab)
#' dim(median_betas)
#' median_betas <- AggregatedEffects(stab, side = "Y")
#' dim(median_betas)
#' }
#'
#' @export
AggregatedEffects <- function(stability, lambda_id = NULL, side = "X", comp = 1,
                              FUN = stats::median, ...) {
  if (!inherits(stability, c("variable_selection", "bi_selection"))) {
    stop("Invalid input for argument 'stability'. This function only applies to the output of VariableSelection() or BiSelection().")
  }

  if (inherits(stability, "variable_selection")) {
    # Extracting corresponding betas over all iterations
    if (ncol(stability$Beta) == stability$params$pk) {
      if (is.null(lambda_id)) {
        if (length(dim(stability$Beta)) == 3) {
          betas <- stability$Beta[ArgmaxId(stability)[1], , ]
        } else {
          betas <- stability$Beta[ArgmaxId(stability)[1], , , ]
        }
      } else {
        if (length(dim(stability$Beta)) == 3) {
          betas <- stability$Beta[lambda_id, , ]
        } else {
          betas <- stability$Beta[lambda_id, , , ]
        }
      }
    } else {
      if (!side %in% c("X", "Y")) {
        warning("Invalid input for argument 'side'. The default value ('X') was used.")
        side <- "X"
      }
      side_id <- grepl(paste0(side, "_"), colnames(stability$Beta))
      comp_id <- grepl(paste0("_PC", comp), colnames(stability$Beta))
      if (is.null(lambda_id)) {
        betas <- stability$Beta[ArgmaxId(stability)[1], side_id & comp_id, ]
      } else {
        betas <- stability$Beta[lambda_id, side_id & comp_id, ]
      }
    }

    # Aggregating the betas conditionally on selection
    if (length(dim(stability$Beta)) == 3) {
      aggregated_betas <- apply(betas, 1, FUN = function(x) {
        FUN(x[x != 0], ...)
      })
      aggregated_betas <- cbind(aggregated_betas)
    } else {
      aggregated_betas <- apply(betas, c(1, 3), FUN = function(x) {
        FUN(x[x != 0], ...)
      })
    }
  } else {
    if (side == "X") {
      betas <- stability$coefX
    } else {
      betas <- stability$coefY
    }
    aggregated_betas <- apply(betas, c(1, 2), FUN = function(x) {
      FUN(x[x != 0], ...)
    })
  }

  # Re-formatting the output
  aggregated_betas[which(is.nan(aggregated_betas))] <- NA

  return(aggregated_betas)
}


#' Calibrated consensus matrix
#'
#' Extracts the (calibrated) consensus matrix. Entries indicate the proportion
#' of subsampling iterations where the two items were in the same cluster, out
#' of all iterations where both items were drawn in the subsample.
#'
#' @param stability output of \code{\link{Clustering}}.
#' @param argmax_id optional matrix of parameter IDs. If \code{argmax_id=NULL},
#'   the calibrated model is used.
#'
#' @return A consensus matrix.
#'
#' @family calibration functions
#'
#' @seealso \code{\link{Clustering}}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of data with clusters
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15),
#'   pk = 10, nu_xc = 1,
#'   ev_xc = 0.95
#' )
#'
#' # Consensus clustering
#' stab <- Clustering(xdata = x, implementation = HierarchicalClustering)
#' ConsensusMatrix(stab)
#' }
#'
#' @export
ConsensusMatrix <- function(stability, argmax_id = NULL) {
  if (class(stability) != "clustering") {
    stop("Invalid input for argument 'stability'. Only applicable to an object of class 'clustering', i.e. to the output of Clustering().")
  }

  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability = stability)
  }
  mat <- stability$coprop[, , argmax_id[1]]

  return(mat)
}


#' Stable cluster membership
#'
#' Generates the (calibrated) stable cluster memberships.
#'
#' @param stability output of \code{\link{Clustering}}.
#' @param linkage character string indicating the type of linkage used in
#'   hierarchical clustering to define the stable clusters. Possible values
#'   include \code{"complete"}, \code{"single"} and \code{"average"} (see
#'   argument \code{"method"} in \code{\link[stats]{hclust}} for a full list).
#' @param argmax_id optional parameter ID. If \code{argmax_id=NULL}, the
#'   calibrated model is used.
#'
#' @return A vector of cluster memberships.
#'
#' @family calibration functions
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(30, 30, 30),
#'   nu_xc = 1
#' )
#' plot(simul)
#'
#' # Consensus clustering
#' stab <- Clustering(
#'   xdata = simul$data
#' )
#' Clusters(stab)
#' }
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


#' Stable attribute weights
#'
#' Creates a boxplots of the distribution of (calibrated) median attribute
#' weights obtained from the COSA algorithm across the subsampling iterations.
#' See examples in \code{\link{Clustering}}.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{Clustering}}.
#' @param at coordinates along the x-axis (more details in
#'   \code{\link[graphics]{boxplot}}).
#' @param col optional vector of colours.
#' @param boxwex box width (more details in \code{\link[graphics]{boxplot}}).
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @param cex.lab font size for labels.
#' @param las orientation of labels on the x-axis (see
#'   \code{\link[graphics]{par}}).
#' @param frame logical indicating if the box around the plot should be drawn
#'   (more details in \code{\link[graphics]{boxplot}}).
#' @param add logical indicating if the boxplot should be added to the current
#'   plot.
#' @param ... additional parameters passed to \code{\link[graphics]{boxplot}}).
#'
#' @return A boxplot.
#'
#' @seealso \code{\link{Clustering}}
#'
#' @export
WeightBoxplot <- function(stability, at = NULL, argmax_id = NULL,
                          col = NULL, boxwex = 0.3,
                          xlab = "", ylab = "Weight", cex.lab = 1.5,
                          las = 3, frame = "F", add = FALSE, ...) {
  # Checking input
  if (all(is.na(stability$Lambda))) {
    stop("Invalid input for argument 'stability'. This function can only be used for weighted clustering.")
  }

  # Defining default colours
  if (is.null(col)) {
    col <- "navy"
  }

  # Extracting ID of calibrated parameters
  if (is.null(argmax_id)) {
    argmax_id <- ArgmaxId(stability)
  }

  # Extracting weights
  y <- stability$Beta[argmax_id[1], , ]

  # Removing zero weights (for sparse methods)
  y[which(y == 0)] <- NA

  xaxt <- "s"
  if (is.null(at)) {
    x <- 1:nrow(y)
  } else {
    x <- at
    xaxt <- "n"
  }

  # Showing the distribution over nonzero weights
  graphics::boxplot(t(y),
    at = x, xlim = range(x),
    col = col, boxcol = col, whiskcol = col,
    staplecol = col,
    whisklty = 1, range = 0, las = las, add = add,
    xlab = xlab, ylab = ylab, cex.lab = cex.lab, frame = frame,
    boxwex = boxwex, xaxt = xaxt, ...
  )
  if (!is.null(at)) {
    graphics::axis(side = 1, at = graphics::axTicks(side = 1))
  }
}


#' Calibration plot
#'
#' Creates a plot showing the stability score as a function of the parameter(s)
#' controlling the level of sparsity in the underlying feature selection
#' algorithm and/or the threshold in selection proportions. See examples in
#' \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#' \code{\link{Clustering}} and \code{\link{BiSelection}}.
#'
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}} or \code{\link{BiSelection}}.
#' @param block_id ID of the block to visualise. Only used for multi-block
#'   stability selection graphical models. If \code{block_id=NULL}, all blocks
#'   are represented in separate panels.
#' @param col vector of colours.
#' @param pch type of point, as in \code{\link[graphics]{points}}.
#' @param cex size of point.
#' @param xlim displayed range along the x-axis. Only used if \code{stability}
#'   is the output of \code{\link{BiSelection}}.
#' @param ylim displayed range along the y-axis. Only used if \code{stability}
#'   is the output of \code{\link{BiSelection}}.
#' @param bty character string indicating if the box around the plot should be
#'   drawn. Possible values include: \code{"o"} (default, the box is drawn), or
#'   \code{"n"} (no box).
#' @param lines logical indicating if the points should be linked by lines. Only
#'   used if \code{stability} is the output of \code{\link{BiSelection}} or
#'   \code{\link{Clustering}}.
#' @param lty line type, as in \code{\link[graphics]{par}}. Only used if
#'   \code{stability} is the output of \code{\link{BiSelection}}.
#' @param lwd line width, as in \code{\link[graphics]{par}}. Only used if
#'   \code{stability} is the output of \code{\link{BiSelection}}.
#' @param show_argmax logical indicating if the calibrated parameter(s) should
#'   be indicated by lines.
#' @param show_pix logical indicating if the calibrated threshold in selection
#'   proportion in \code{X} should be written for each point. Only used if
#'   \code{stability} is the output of \code{\link{BiSelection}}.
#' @param show_piy logical indicating if the calibrated threshold in selection
#'   proportion in \code{Y} should be written for each point. Only used if
#'   \code{stability} is the output of \code{\link{BiSelection}} with
#'   penalisation of the outcomes.
#' @param offset distance between the point and the text, as in
#'   \code{\link[graphics]{text}}. Only used if \code{show_pix=TRUE} or
#'   \code{show_piy=TRUE}.
#' @param legend logical indicating if the legend should be included.
#' @param legend_length length of the colour bar. Only used if \code{stability}
#'   is the output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}.
#' @param legend_range range of the colour bar. Only used if \code{stability} is
#'   the output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}.
#' @param ncol integer indicating the number of columns in the legend.
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @param zlab label of the z-axis. Only used if \code{stability} is the output
#'   of \code{\link{VariableSelection}} or \code{\link{GraphicalModel}}.
#' @param xlas orientation of labels on the x-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylas orientation of labels on the y-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param zlas orientation of labels on the z-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param cex.lab font size for labels.
#' @param cex.axis font size for axes.
#' @param cex.legend font size for text legend entries.
#' @param xgrid logical indicating if a vertical grid should be drawn. Only used
#'   if \code{stability} is the output of \code{\link{BiSelection}}.
#' @param ygrid logical indicating if a horizontal grid should be drawn. Only
#'   used if \code{stability} is the output of \code{\link{BiSelection}}.
#' @param params vector of possible parameters if \code{stability} is of class
#'   \code{bi_selection}. The order of these parameters defines the order in
#'   which they are represented. Only used if \code{stability} is the output of
#'   \code{\link{BiSelection}}.
#'
#' @return A calibration plot.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{Clustering}}, \code{\link{BiSelection}}
#'
#' @export
CalibrationPlot <- function(stability, block_id = NULL,
                            col = NULL,
                            pch = 19, cex = 0.7,
                            xlim = NULL, ylim = NULL, bty = "o",
                            lines = TRUE, lty = 3, lwd = 2,
                            show_argmax = TRUE,
                            show_pix = FALSE, show_piy = FALSE, offset = 0.3,
                            legend = TRUE, legend_length = NULL, legend_range = NULL, ncol = 1,
                            xlab = NULL, ylab = NULL, zlab = expression(italic(q)),
                            xlas = 2, ylas = NULL, zlas = 2,
                            cex.lab = 1.5, cex.axis = 1, cex.legend = 1.2,
                            xgrid = FALSE, ygrid = FALSE,
                            params = c("ny", "alphay", "nx", "alphax")) {
  oldpar <- graphics::par("xpd", "xaxs", "yaxs", no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # To deal with later: showing calibration of clustering or selection
  # clustering <- FALSE
  clustering <- ifelse(inherits(stability, "clustering"), yes = TRUE, no = FALSE)
  heatmap <- TRUE

  if (clustering) {
    ylas <- 1
    CalibrationCurve(
      stability = stability, bty = bty, xlab = xlab, ylab = ylab,
      cex.lab = cex.lab, cex.axis = cex.axis, cex.legend = cex.legend,
      pch = pch, lines = lines, col = col, legend = legend, ncol = ncol
    )
  } else {
    ylas <- 0

    if (inherits(stability, "bi_selection")) {
      # Extracting summary information
      x <- stability$summary_full

      # Checking input
      params <- unique(params)
      all_params <- colnames(stability$summary)
      all_params <- all_params[!all_params %in% c("comp", "S", "pix", "piy")]
      if (any(!all_params %in% params)) {
        params <- all_params
        warning(paste0(
          "Invalid input for argument 'params'. Please provide a vector with all the following: ",
          paste(all_params, collapse = ", "), "."
        ))
      }
      params <- params[params %in% all_params]

      # Identifying parameters
      params <- params[params %in% colnames(x)]

      # Defining default arguments
      if (is.null(ylab)) {
        ylab <- "Stability Score"
      }

      if (is.null(xlab)) {
        if (length(params) > 1) {
          xlab <- ""
        } else {
          xlab <- expression(n[X])
        }
      }

      if (is.null(col)) {
        col <- grDevices::colorRampPalette(c("navy", "darkred"))(nrow(stability$summary))
      } else {
        col <- grDevices::colorRampPalette(col)(nrow(stability$summary))
      }

      if (length(unique(x$comp)) == 1) {
        legend <- FALSE
      }

      if (is.null(xlim)) {
        xlim <- c(0.5, max(sapply(split(x, f = x$comp), nrow)) + 0.5)
      }

      if (is.null(ylim)) {
        ylim <- range(x$S)
        if (legend) {
          ylim[2] <- ylim[2] + diff(ylim) * 0.15
        }
      }

      # Drawing one set of points per component
      for (comp_id in unique(x$comp)) {
        tmp <- x[which(x$comp == comp_id), ]

        # Ensuring increasing ny
        tmp <- tmp[do.call(order, tmp[, params, drop = FALSE]), ]
        # if ("ny" %in% colnames(tmp)) {
        #   tmp=tmp[order(tmp$ny, tmp$nx), ]
        # }
        # tmp=tmp[order(lapply(params, FUN=function(param_id){with(tmp, eval(parse(text=param_id)))})),]
        # if ("alphax" %in% colnames(tmp))

        if (comp_id == min(x$comp)) {
          # Initialising the plot
          plot(NA,
            xlim = xlim, ylim = ylim, bty = bty,
            xlab = xlab, ylab = ylab, cex.lab = cex.lab,
            cex.axis = cex.axis,
            xaxt = "n", las = ylas
          )

          # Defining vertical grid
          if (xgrid) {
            withr::local_par(list(xpd = FALSE))
            graphics::abline(v = 1:nrow(tmp), lty = 3, col = "grey")
          }

          # Defining horizontal grid
          if (ygrid) {
            withr::local_par(list(xpd = FALSE))
            graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
          }

          # Adding x-axes
          for (param_id in 1:length(params)) {
            if (param_id == 1) {
              graphics::axis(
                side = 1, at = 1:nrow(tmp),
                labels = tmp[, rev(params)[param_id]],
                cex.axis = cex.axis, las = xlas
              )
            } else {
              ids <- c(1, which(diff(tmp[, rev(params)[param_id]]) != 0) + 1)
              ids <- c(ids - 0.5, nrow(tmp) + 0.5)
              graphics::axis(side = 1, at = ids, labels = NA, line = (param_id - 1) * 3)
              withr::local_par(list(xpd = FALSE))
              graphics::abline(v = ids, lty = 2)
              ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
              graphics::axis(
                side = 1, at = ids, labels = tmp[ids, rev(params)[param_id]],
                line = (param_id - 1) * 3, tick = FALSE, cex.axis = cex.axis, las = xlas
              )
            }
          }
          # graphics::axis(side = 1, at = 1:nrow(tmp), labels = tmp$nx, cex.axis = cex.axis, las = xlas)
          # if ("ny" %in% colnames(tmp)) {
          #   ids <- c(which(!duplicated(tmp$ny)) - 0.5, nrow(tmp) + 0.5)
          #   graphics::axis(side = 1, at = ids, labels = NA, line = 3)
          #   withr::local_par(list(xpd = FALSE))
          #   graphics::abline(v = ids, lty = 2)
          #   ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
          #   graphics::axis(side = 1, at = ids, labels = unique(tmp$ny), line = 3, tick = FALSE, cex.axis = cex.axis, las = xlas)
          # }

          # Adding x-labels
          if (length(params) > 1) {
            for (param_id in 1:length(params)) {
              if (rev(params)[param_id] == "nx") {
                graphics::mtext(text = expression(n[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
              }
              if (rev(params)[param_id] == "alphax") {
                graphics::mtext(text = expression(alpha[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
              }
              if (rev(params)[param_id] == "ny") {
                graphics::mtext(text = expression(n[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
              }
              if (rev(params)[param_id] == "alphay") {
                graphics::mtext(text = expression(alpha[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
              }
            }
            # graphics::mtext(text = expression(n[X]), side = 1, at = -nrow(tmp) * 0.06, line = 1, cex = cex.lab)
            # graphics::mtext(text = expression(n[Y]), side = 1, at = -nrow(tmp) * 0.06, line = 4, cex = cex.lab)
          }
        }

        # Adding calibrated lines
        if (show_argmax) {
          withr::local_par(list(xpd = FALSE))
          graphics::abline(v = which.max(tmp$S), lty = 3, col = col[comp_id])
        }

        # Adding lines
        if (lines) {
          # if ("ny" %in% colnames(tmp)) {
          #   for (y_value in unique(tmp$ny)) {
          #     graphics::lines(which(tmp$ny == y_value),
          #                     tmp[which(tmp$ny == y_value), "S"],
          #                     col = col[comp_id],
          #                     lty = lty, lwd = lwd
          #     )
          #   }
          # } else {
          graphics::lines(1:nrow(tmp),
            tmp$S,
            col = col[comp_id],
            lty = lty, lwd = lwd
          )
          # }
        }

        # Adding data points
        graphics::points(tmp$S,
          pch = pch,
          col = col[comp_id],
          cex = cex
        )

        # Adding pi values
        if ((show_pix) & (!show_piy)) {
          graphics::text(tmp$S,
            labels = tmp$pix,
            col = col[comp_id],
            cex = cex, pos = 3,
            offset = offset
          )
        }

        if ((!show_pix) & (show_piy)) {
          graphics::text(tmp$S,
            labels = tmp$piy,
            col = col[comp_id],
            cex = cex, pos = 3,
            offset = offset
          )
        }

        if ((show_pix) & (show_piy)) {
          for (k in 1:nrow(tmp)) {
            graphics::text(k, tmp[k, "S"],
              labels = eval(parse(text = paste0("expression(pi[x]*' = ", tmp[k, "pix"], " ; '*pi[y]*' = ", tmp[k, "piy"], "')"))),
              col = col[comp_id],
              cex = cex, pos = 3,
              offset = offset
            )
          }
        }
      }

      # Adding legend
      if (legend) {
        graphics::legend("top",
          col = col, lty = lty, pch = pch, lwd = lwd,
          legend = paste0("Component ", unique(x$comp)),
          horiz = TRUE, bg = "white", cex = cex.legend
        )
      }
    } else {
      # Defining default arguments
      if (heatmap) {
        metric <- "both"
        if (is.null(col)) {
          col <- c("ivory", "navajowhite", "tomato", "darkred")
        }
        if (is.null(ylab)) {
          if (clustering) {
            ylab <- expression(n[c])
          } else {
            ylab <- expression(pi)
          }
        }
      } else {
        metric <- "lambda"
        if (is.null(col)) {
          col <- "navy"
        }
        if (is.null(ylab)) {
          ylab <- "Stability Score"
        }
      }
      if (is.null(xlab)) {
        xlab <- expression(lambda)
      }

      # Extracting the number of blocks/components
      if ((stability$methods$type == "graphical_model") & (is.null(block_id))) {
        bigblocks <- BlockMatrix(stability$params$pk)
        bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
        N_blocks <- unname(table(bigblocks_vect))
        blocks <- unique(as.vector(bigblocks_vect))
        names(N_blocks) <- blocks
        nblocks <- max(blocks)
        block_id <- 1:nblocks
      } else {
        block_id <- 1
      }
      nblocks <- length(block_id)

      if (metric == "both") {
        for (b in block_id) {
          # Extracting the stability scores
          if (clustering) {
            mat <- matrix(stability$Sc, ncol = length(unique(stability$Lambda)))
            rownames(mat) <- formatC(unique(stability$nc), format = "f", digits = 0)
            colnames(mat) <- formatC(unique(stability$Lambda), format = "e", digits = 2)
            mat <- t(mat)
            mat <- mat[nrow(mat):1, ncol(mat):1]
            ids <- which(apply(mat, 1, FUN = function(x) {
              any(!is.na(x))
            }))
            mat <- mat[ids, , drop = FALSE]
          } else {
            if (length(stability$params$pk) == 1) {
              mat <- stability$S_2d
              ids <- which(apply(mat, 1, FUN = function(x) {
                any(!is.na(x))
              }))
              mat <- mat[ids, , drop = FALSE]
            } else {
              mat <- stability$S_2d[, , b]
              ids <- which(apply(mat, 1, FUN = function(x) {
                any(!is.na(x))
              }))
              mat <- mat[ids, , drop = FALSE]
            }

            # Setting row and column names
            colnames(mat) <- stability$params$pi_list
            if (grepl("penalised", tolower(stability$methods$implementation))) {
              rownames(mat) <- formatC(stability$Lambda[, b], format = "e", digits = 2)[ids]
            } else {
              rownames(mat) <- (stability$Lambda[, b])[ids]
            }
          }

          # Extracting corresponding numbers of selected variables (q)
          Q <- stability$Q[, b]
          Q <- Q[ids]

          # Heatmap representation
          Heatmap(t(mat[nrow(mat):1, ncol(mat):1]),
            col = col, bty = bty, axes = FALSE,
            text = FALSE, cex = 1, digits = 2,
            legend = legend, legend_length = legend_length, legend_range = legend_range
          )

          # Adding calibrated lines
          if (show_argmax) {
            withr::local_par(list(xpd = FALSE))
            if (clustering) {
              myargmax <- Argmax(stability)
              print(myargmax)
              argmax_id <- c(
                which(rownames(t(mat[nrow(mat):1, ncol(mat):1])) == formatC(myargmax[1], format = "f", digits = 0)),
                which(colnames(t(mat[nrow(mat):1, ncol(mat):1])) == formatC(myargmax[2], format = "e", digits = 2))
              )
              print(argmax_id)
              graphics::abline(h = ncol(mat) - argmax_id[1] + 0.5, lty = 3)
              graphics::abline(v = argmax_id[2] - 0.5, lty = 3)
            } else {
              graphics::abline(v = nrow(mat) - which(stability$Lambda[ids, b] == Argmax(stability)[b, 1]) + 0.5, lty = 3)
              graphics::abline(h = which.min(abs(as.numeric(colnames(mat)) - Argmax(stability)[b, 2])) - 0.5, lty = 3)
            }
          }

          # Including axes
          if (clustering) {
            graphics::axis(
              side = 2, at = (1:ncol(mat)) - 0.5, las = ylas, cex.axis = cex.axis,
              labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 0)
            )
          } else {
            graphics::axis(
              side = 2, at = (1:ncol(mat)) - 0.5, las = ylas, cex.axis = cex.axis,
              labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2)
            )
          }
          if (grepl("penalised", tolower(stability$methods$implementation))) {
            graphics::axis(
              side = 3, at = (1:nrow(mat)) - 0.5, las = zlas, cex.axis = cex.axis,
              labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0))
            )
            graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = xlas, labels = rev(rownames(mat)), cex.axis = cex.axis)
          } else {
            graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = xlas, labels = rev(rownames(mat)), cex.axis = cex.axis)
          }

          # Including axis labels
          graphics::mtext(text = ylab, side = 2, line = 3.5, cex = cex.lab)
          graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
          if (grepl("penalised", tolower(stability$methods$implementation))) {
            graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
          }
        }
      } else {
        if (metric == "lambda") {
          for (b in block_id) {
            # Extracting the stability scores
            if (length(stability$params$pk) == 1) {
              mat <- stability$S_2d
              ids <- which(apply(mat, 1, FUN = function(x) {
                any(!is.na(x))
              }))
              mat <- mat[ids, , drop = FALSE]
            } else {
              mat <- stability$S_2d[, , b]
              ids <- which(apply(mat, 1, FUN = function(x) {
                any(!is.na(x))
              }))
              mat <- mat[ids, , drop = FALSE]
            }

            # Extracting the best stability score (with optimal pi) for each lambda value
            vect <- apply(mat, 1, max, na.rm = TRUE)

            # Extracting corresponding numbers of selected variables (q)
            Q <- stability$Q[, b, drop = FALSE]
            Q <- Q[ids]

            # Extracting corresponding lambda values
            Lambda <- stability$Lambda[ids, b, drop = FALSE]

            # Re-ordering by decreasing lambda
            ids <- sort.list(Lambda, decreasing = TRUE)
            Lambda <- Lambda[ids]
            Q <- Q[ids]
            vect <- vect[ids]

            if (is.null(xlim)) {
              xlim <- range(Lambda, na.rm = TRUE)
            }

            if (is.null(ylim)) {
              ylim <- range(vect)
            }

            # Initialising the plot
            plot(NA,
              xlim = xlim, ylim = ylim, bty = bty,
              xlab = "", ylab = ylab, cex.lab = cex.lab,
              cex.axis = cex.axis,
              xaxt = "n", las = ylas
            )

            # Defining horizontal grid
            if (ygrid) {
              withr::local_par(list(xpd = FALSE))
              graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
            }

            # Adding calibrated lines
            if (show_argmax) {
              withr::local_par(list(xpd = FALSE))
              graphics::abline(h = max(vect), lty = 3, col = col[1])
              graphics::abline(v = Lambda[which.max(vect)], lty = 3, col = col[1])
            }

            # Adding lines
            if (lines) {
              graphics::lines(Lambda, vect, col = col[1], lty = lty, lwd = lwd)
            }

            # Adding data points
            graphics::points(Lambda, vect, pch = pch, col = col[1], cex = cex)

            # Adding x-axis and z-axis and their labels
            lseq <- grDevices::axisTicks(range(Lambda, na.rm = TRUE), log = FALSE)
            xseq <- 1
            for (i in 1:length(lseq)) {
              xseq <- c(xseq, which.min(abs(Lambda - lseq[i])))
            }
            xseq <- c(xseq, length(Lambda))
            xseq <- unique(xseq)
            if (xgrid) {
              withr::local_par(list(xpd = FALSE))
              graphics::abline(v = Lambda[xseq], lty = 3, col = "grey")
            }
            graphics::axis(side = 1, at = Lambda[xseq], labels = formatC(Lambda[xseq], format = "e", digits = 2), las = xlas, cex.axis = cex.axis)
            graphics::axis(
              side = 3, at = Lambda[xseq], las = xlas,
              labels = formatC(Q[xseq], format = "f", big.mark = ",", digits = 0), cex.axis = cex.axis
            )
            graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
            graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
          }
        }
      }
    }
  }
}


#' Calibration curve (internal)
#'
#' Creates a calibration curve for consensus clustering.
#'
#' @inheritParams CalibrationPlot
#'
#' @return a calibration curve.
#'
#' @keywords internal
CalibrationCurve <- function(stability,
                             bty = "o",
                             xlab = NULL,
                             ylab = NULL,
                             cex.axis = 1,
                             cex.lab = 1.5,
                             cex.legend = 1.2,
                             pch = 19,
                             lines = TRUE,
                             col = NULL,
                             legend = TRUE,
                             ncol = 1) {
  # Checking inputs
  if (is.null(col)) {
    col <- c("navy", "forestgreen", "tomato")
  }
  if (is.null(xlab)) {
    xlab <- expression(n[C])
  }
  if (is.null(ylab)) {
    ylab <- "Consensus score"
  }

  # Creating objects
  y <- stability$Sc
  x <- stability$nc
  if (any(!is.na(stability$Lambda))) {
    z <- round(stability$Lambda, digits = 5)
  } else {
    z <- rep(0, length(stability$nc))
  }

  # Preparing colours
  if (length(unique(z)) > 1) {
    mycolours <- grDevices::colorRampPalette(col)(length(unique(z)))
  } else {
    mycolours <- col[1]
  }
  names(mycolours) <- unique(z)

  # Initialising plot
  plot(NA,
    xlim = c(0, max(stability$nc)), ylim = c(0, 1),
    xlab = xlab, ylab = ylab, cex.axis = cex.axis,
    las = 1, cex.lab = cex.lab, bty = bty
  )

  # Adding lines
  for (lambda in unique(z)) {
    ids <- which(z == lambda)
    graphics::points(x[ids], y[ids], pch = 18, col = mycolours[as.character(lambda)])
    if (lines) {
      graphics::lines(x[ids], y[ids], lty = 1, lwd = 0.5, col = mycolours[as.character(lambda)])
    }
  }

  # Identifying calibrated parameter(s)
  graphics::abline(v = Argmax(stability)[1], lty = 2, col = "darkred")

  # Adding legend
  if (any(!is.na(stability$Lambda))) {
    if (legend) {
      if (length(unique(stats::na.exclude(stability$Q))) == 1) {
        legend("topright",
          legend = unique(formatC(stability$Lambda, format = "f", digits = 2)),
          pch = 15, col = mycolours, cex = cex.legend, bty = "n", title = expression(lambda), ncol = ncol
        )
      } else {
        legend("topright",
          legend = paste0(
            unique(formatC(stability$Lambda[, 1], format = "f", digits = 2)),
            " (", unique(stability$Q[, 1]), ")"
          ),
          pch = 15, col = mycolours, cex = cex.legend, bty = "n", title = expression(lambda), ncol = ncol
        )
      }
    }
  }
}

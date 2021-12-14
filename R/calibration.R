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
#'   provided. This argument cannot be used with \code{clustering=TRUE}.
#' @param clustering logical indicating whether indices of calibrated clustering
#'   parameters (\code{clustering=TRUE}) or variable selection
#'   (\code{clustering=FALSE}) should be extracted. This argument is only used
#'   if \code{stability} is the output from \code{\link{Clustering}}.
#'
#' @return A matrix of parameter indices. In multi-block graphical modelling,
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
ArgmaxId <- function(stability = NULL, S = NULL, clustering = FALSE) {
  if ((is.null(stability)) & (is.null(S))) {
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }
  if (clustering) {
    if (!is.null(S)) {
      stop("Invalid input. Argument 'stability' needs to be supplied with clustering = TRUE.")
    }
  }
  if (is.null(S)) {
    if (clustering) {
      if (length(unique(stability$Lambda)) > 1) {
        # Identifying best number of contributing variables
        lambda_hat <- stability$Lambda[which.max(stability$S), 1]
        ids <- which(as.character(stability$Lambda) == lambda_hat)
      } else {
        ids <- 1:nrow(stability$Sc)
      }
      Sc <- stability$Sc[ids, 1]
      Sc_2d <- stability$Sc_2d[ids, , drop = FALSE]

      # Identifying best number of clusters
      argmax_id <- matrix(NA, nrow = 1, ncol = 2)
      id <- which.max(Sc)
      argmax_id[, 1] <- ids[id]
      tmpSc <- Sc_2d[id, ]
      argmax_id[, 2] <- which.max(tmpSc)
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
  colnames(argmax_id) <- c("lambda_id", "pi_id")
  return(argmax_id)
}


#' Calibrated parameters
#'
#' Extracts calibrated parameter values in stability selection.
#'
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}.
#' @inheritParams ArgmaxId
#'
#' @return A matrix of parameter values. The first column (\code{lambda}])
#'   denotes the calibrated hyper-parameter of the underlying algorithm. The
#'   second column (\code{pi}) is the calibrated threshold in
#'   selection/co-membership proportions. In multi-block graphical modelling,
#'   rows correspond to different blocks.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
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
#' args <- Argmax(stab)
#'
#'
#' ## Clustering
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15)
#' )
#'
#' # Consensus clustering
#' stab <- Clustering(xdata = simul$data, Lambda = c(1.1, 1.5, 2))
#'
#' # Extracting calibrated parameters for variable selection
#' Argmax(stab, clustering = FALSE)
#'
#' # Extracting calibrated parameters for co-membership
#' Argmax(stab, clustering = TRUE)
#' }
#'
#' @export
Argmax <- function(stability, clustering = FALSE) {
  argmax <- matrix(NA, nrow = ncol(stability$Lambda), ncol = 2)
  if (clustering) {
    id <- ArgmaxId(stability = stability, clustering = clustering)
    argmax[, 1] <- stability$nc[id[1], 1]
    argmax[, 2] <- stability$params$pi_list[id[2]]
  } else {
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
  if (class(stability) == "bi_selection") {
    if ("selectedY" %in% names(stability)) {
      A <- Square(cbind(stability$selectedX, stability$selectedY))
    } else {
      A <- Square(stability$selectedX)
    }
  } else {
    if (class(stability) == "graphical_model") {
      A <- matrix(0, ncol = ncol(stability$selprop), nrow = nrow(stability$selprop))
    } else {
      A <- matrix(0, ncol = ncol(stability$coprop), nrow = nrow(stability$coprop))
    }
    bigblocks <- BlockMatrix(stability$params$pk)
    if (is.null(argmax_id)) {
      if (class(stability) == "graphical_model") {
        argmax_id <- ArgmaxId(stability)
        argmax <- Argmax(stability)
      } else {
        argmax_id <- ArgmaxId(stability, clustering = TRUE)
        argmax <- Argmax(stability, clustering = TRUE)
      }
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
      if (class(stability) == "graphical_model") {
        A_block <- ifelse(stability$selprop[, , argmax_id[block_id, 1]] >= argmax[block_id, 2], 1, 0)
      } else {
        A_block <- ifelse(stability$coprop[, , argmax_id[block_id, 1]] >= argmax[block_id, 2], 1, 0)
      }
      A_block[lower.tri(A_block)] <- 0
      A_block <- A_block + t(A_block) # for symmetry
      if (length(stability$params$pk) > 1) {
        A_block[bigblocks != block_id] <- 0
      }
      A <- A + A_block
    }
  }
  A[is.na(A)] <- 0
  return(A)
}


#' Set of stably selected variables
#'
#' Builds the (calibrated) set of stably selected variables.
#'
#' @inheritParams Adjacency
#' @param stability output of \code{\link{VariableSelection}},
#'   or \code{\link{BiSelection}}.
#'
#' @return A binary vector encoding the selection status of the variables
#'   (\code{1} if selected, \code{0} otherwise).
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{BiSelection}}
#'
#' @examples
#' \dontrun{
#'
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
#' myids <- matrix(c(80, 10), nrow = 1)
#' stab$Lambda[myids[1], 1] # corresponding penalty
#' stab$params$pi_list[myids[2]] # corresponding threshold
#' selected <- SelectedVariables(stab, argmax_id = myids)
#' }
#' @export
SelectedVariables <- function(stability, argmax_id = NULL) {
  if (!class(stability) %in% c("clustering", "variable_selection", "bi_selection")) {
    stop("Invalid input for argument 'stability'. This function only applies to outputs from VariableSelection() or BiSelection().")
  }

  if (class(stability) == "clustering") {
    selprop <- SelectionProportions(stability, argmax_id = argmax_id)
    if (any(selprop != 1)) {
      score <- StabilityScore(
        selprop = selprop,
        K = stability$params$K,
        pi_list = stability$params$pi_list,
        n_cat = stability$params$n_cat
      )
      stability_selected <- ifelse(selprop >= stability$params$pi_list[which.max(score)],
        yes = 1, no = 0
      )
    } else {
      stability_selected <- selprop
    }
  }

  if (class(stability) == "variable_selection") {
    if (is.null(argmax_id)) {
      argmax_id <- ArgmaxId(stability)
    }
    stability_selected <- ifelse(stability$selprop[argmax_id[1], ] >= stability$params$pi_list[argmax_id[2]],
      yes = 1, no = 0
    )
  }

  if (class(stability) == "bi_selection") {
    stability_selected <- stability$selectedX
  }

  return(stability_selected)
}


#' Stable cluster membership
#'
#' Builds the (calibrated) stable clusters as connected components of the graph
#' defined from stable co-membership.
#'
#' @inheritParams Adjacency
#' @param adjacency adjacency matrix or output of \code{\link{GraphicalModel}}
#'   or \code{\link{Clustering}}.
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
#' # Consensus clustering
#' stab <- GraphicalModel(
#'   xdata = simul$data,
#'   Lambda = seq(2, ncol(simul$data)),
#'   implementation = HierarchicalClustering
#' )
#'
#' # Stable cluster membership
#' groups <- Clusters(stab)
#'
#' # Network representation of stable co-membership
#' set.seed(1)
#' plot(Graph(CoMembership(groups),
#'   satellites = TRUE,
#'   node_colour = distinctColorPalette(length(unique(groups)))[groups]
#' ))
#' }
#' @export
Clusters <- function(adjacency = NULL, argmax_id = NULL) {
  if (!is.matrix(adjacency)) {
    # Computing stable co-membership matrix
    adjacency <- Adjacency(stability = adjacency, argmax_id = argmax_id)
  }

  # Extracting stable connected components
  mymembership <- igraph::components(Graph(adjacency, satellites = TRUE))$membership

  return(mymembership)
}


#' Selection proportions
#'
#' Extracts the selection proportions of the (calibrated) stability selection
#' model or co-membership proportions of the (calibrated) consensus clustering
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
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
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
#'
#' # Sparse PLS
#' stab <- BiSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
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

  if (class(stability) == "graphical_model") {
    out <- SelectionProportionsGraphical(stability = stability, argmax_id = argmax_id)
  }
  if (class(stability) == "variable_selection") {
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (class(stability) == "clustering") {
    argmax_id <- ArgmaxId(stability)
    out <- SelectionProportionsRegression(stability = stability, argmax_id = argmax_id)
  }
  if (class(stability) == "bi_selection") {
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
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#'
#' # Coefficients of visited models
#' coefs <- focus:::Coefficients(stab)
#' dim(coefs)
#'
#' # Coefficients of the first fitted model
#' coefs <- focus:::Coefficients(stab, iterations = 1)
#' dim(coefs)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   implementation = SparsePLS, family = "gaussian"
#' )
#'
#' # Coefficients of visited models
#' coefs <- focus:::Coefficients(stab, side = "Y", )
#' dim(coefs)
#' }
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


#' Average coefficients conditionally on selection
#'
#' Extracts the average coefficients of the (calibrated) models conditionally on
#' selection across resampling iterations. This function can be applied to the
#' output of \code{\link{VariableSelection}}.
#'
#' @param stability output of \code{\link{VariableSelection}}.
#' @param side character string indicating if coefficients of the predictor
#'   (\code{side="X"}) or outcome (\code{side="Y"}) coefficients should be
#'   returned. Option \code{side="Y"} is only applicable to PLS models.
#' @param comp component ID. Only applicable to PLS models.
#' @param lambda_id parameter ID with respect to the grid \code{Lambda}. If
#'   \code{NULL}, average coefficients across the models run with the calibrated
#'   parameter are returned.
#' @param FUN function to use to aggregate coefficients of visited models over
#'   resampling iterations. Recommended functions include \code{\link[stats]{median}}
#'   or \code{\link[base]{mean}}.
#' @param ... additional arguments to be passed to \code{FUN}.
#'
#' @return A matrix of average coefficients across resampling iterations.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Recalibrate}}
#'
#' @examples
#' \dontrun{
#'
#' # Example with univariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' median_betas <- AggregatedEffects(stab)
#'
#' # Comparison with recalibrated model
#' recalibrated <- Recalibrate(xdata = simul$xdata, ydata = simul$ydata, stability = stab)
#' recalib_betas <- recalibrated$coefficients[-1]
#' plot(median_betas[names(recalib_betas), ], recalib_betas,
#'   panel.first = abline(0, 1, lty = 2)
#' )
#'
#' # Regression with multivariate outcomes
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = c(20, 30), family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "mgaussian")
#' median_betas <- AggregatedEffects(stab)
#' dim(median_betas)
#'
#' # Sparse PLS with multivariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
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
#' @export
AggregatedEffects <- function(stability, lambda_id = NULL, side = "X", comp = 1,
                             FUN = stats::median, ...) {
  if (!class(stability) %in% c("variable_selection", "bi_selection")) {
    stop("Invalid input for argument 'stability'. This function only applies to the output of VariableSelection() or BiSelection().")
  }

  if (class(stability) == "variable_selection") {
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
      av_betas <- apply(betas, 1, FUN = function(x) {
        FUN(x[x != 0], ...)
      })
      av_betas <- cbind(av_betas)
    } else {
      av_betas <- apply(betas, c(1, 3), FUN = function(x) {
        FUN(x[x != 0], ...)
      })
    }
  } else {
    if (side == "X") {
      betas <- stability$coefX
    } else {
      betas <- stability$coefY
    }
    av_betas <- apply(betas, c(1, 2), FUN = function(x) {
      FUN(x[x != 0], ...)
    })
  }

  # Re-formatting the output
  av_betas[which(is.nan(av_betas))] <- NA

  return(av_betas)
}

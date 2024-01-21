#' Transforms NA into NULL
#'
#' Returns a vector with no missing values or NULL if there are no non-missing
#' values.
#'
#' @param x input vector.
#'
#' @return A vector without missing values or NULL.
#'
#' @keywords internal
NAToNULL <- function(x) {
  if (any(!is.na(x))) {
    return(x = x[!is.na(x)])
  } else {
    return(NULL)
  }
}


#' Adjacency matrix from object
#'
#' Returns the adjacency matrix from an
#' \code{\link[igraph:igraph-package]{igraph}} object or from the output of
#' simulation or inference functions from the present package.
#'
#' @param object input object.
#'
#' @return A vector without missing values or NULL.
#'
#' @keywords internal
AdjacencyFromObject <- function(object) {
  if (inherits(object, "matrix")) {
    theta <- object
  } else {
    # From igraph object
    if (inherits(object, "igraph")) {
      theta <- as.matrix(igraph::get.adjacency(object))
    }

    # From output of VariableSelection()
    if (inherits(object, "variable_selection")) {
      theta <- cbind(SelectedVariables(object))
      if (ncol(theta) == 1) {
        colnames(theta) <- "outcome"
      }
      theta <- Square(theta)
    }

    # From output of GraphicalModel()
    if (inherits(object, "graphical_model")) {
      theta <- Adjacency(object)
    }

    # From output of BiSelection()
    if (inherits(object, "bi_selection")) {
      if ("selected" %in% names(object)) {
        selected <- object$selected # PLS
      } else {
        selected <- object$selectedX # PCA
      }
      theta <- Square(selected)
    }

    # From output of SimulateRegression() or SimulateComponents()
    if (inherits(object, c("simulation_regression", "simulation_components"))) {
      theta <- cbind(object$theta)
      if (ncol(theta) == 1) {
        colnames(theta) <- "outcome"
      }
      theta <- Square(theta)
    }

    # From output of SimulateGraphical()
    if (inherits(object, "simulation_graphical_model")) {
      theta <- object$theta
    }
  }

  return(theta)
}


#' Adjacency from bipartite
#'
#' Generates a symmetric adjacency matrix encoding a bipartite graph.
#'
#' @param x matrix encoding the edges between two types of nodes (rows and
#'   columns).
#'
#' @return A symmetric adjacency matrix encoding a bipartite graph.
#'
#' @examples
#' # Simulated links between two sets
#' set.seed(1)
#' mat <- matrix(sample(c(0, 1), size = 15, replace = TRUE),
#'   nrow = 5, ncol = 3
#' )
#'
#' # Adjacency matrix of a bipartite graph
#' Square(mat)
#' @export
Square <- function(x) {
  # Assigning row and column names
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("row", seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("col", seq_len(ncol(x)))
  }

  # Defining the square matrix
  adjacency <- rbind(
    cbind(matrix(0, nrow = nrow(x), ncol = nrow(x)), x),
    cbind(t(x), matrix(0, nrow = ncol(x), ncol = ncol(x)))
  )
  rownames(adjacency) <- colnames(adjacency) <- c(rownames(x), colnames(x))

  return(adjacency)
}


#' Categorical from dummy variables
#'
#' Generates a single categorical variable from corresponding dummy variables.
#'
#' @param x matrix of dummy variables.
#' @param verbose logical indicating if messages should be printed.
#'
#' @return A single categorical variable (numeric).
#'
#' @keywords internal
DummyToCategories <- function(x, verbose = FALSE) {
  x_original <- x
  x <- matrix(0, nrow = nrow(x_original), ncol = ncol(x_original))
  for (j in seq_len(ncol(x))) {
    tmp <- as.factor(x_original[, j])
    if (verbose) {
      message(paste0("Reference category for column ", j, ": ", levels(tmp)[1]))
      message(paste0("Other category for column ", j, ": ", levels(tmp)[2]))
    }
    x[, j] <- (as.numeric(tmp) - 1) * j
  }
  x <- apply(x, 1, sum)
  return(x)
}


#' Pairwise co-membership
#'
#' Generates a symmetric and binary matrix indicating, if two items are
#' co-members, i.e. belong to the same cluster.
#'
#' @param groups vector of group membership.
#'
#' @return A symmetric and binary matrix.
#'
#' @examples
#' # Simulated grouping structure
#' mygroups <- c(rep(1, 3), rep(2, 5), rep(3, 2))
#'
#' # Co-membership matrix
#' CoMembership(mygroups)
#' @export
CoMembership <- function(groups) {
  if (length(unique(groups)) > 1) {
    # Building binary cluster membership for each feature
    V <- stats::model.matrix(~ as.factor(groups) - 1)

    # Building cluster co-membership
    comembership <- V %*% t(V)
  } else {
    comembership <- matrix(1, nrow = length(groups), ncol = length(groups))
  }

  # Re-formatting co-membership matrix
  diag(comembership) <- 1
  rownames(comembership) <- colnames(comembership) <- names(groups)

  return(comembership)
}


#' Concatenate stability objects
#'
#' Generates a single stability object from two stability objects. This function
#' is used to concatenate results when using \code{\link[nloptr]{nloptr}}.
#'
#' @param stability1 a stability object.
#' @param stability2 another stability object (optional).
#'
#' @return A single stability object.
#'
#' @keywords internal
Concatenate <- function(stability1, stability2 = NULL, order_output = FALSE) {
  if (!is.null(stability2)) {
    stability1$S <- rbind(stability1$S, stability2$S)
    stability1$Lambda <- rbind(stability1$Lambda, stability2$Lambda)
    stability1$Q <- rbind(stability1$Q, stability2$Q)
    stability1$Q_s <- rbind(stability1$Q_s, stability2$Q_s)
    stability1$P <- rbind(stability1$P, stability2$P)
    stability1$PFER <- rbind(stability1$PFER, stability2$PFER)
    stability1$FDP <- rbind(stability1$FDP, stability2$FDP)
    stability1$S_2d <- rbind(stability1$S_2d, stability2$S_2d)
    stability1$PFER_2d <- rbind(stability1$PFER_2d, stability2$PFER_2d)
    stability1$FDP_2d <- rbind(stability1$FDP_2d, stability2$FDP_2d)
    if (stability1$methods$type == "variable_selection") {
      stability1$selprop <- rbind(stability1$selprop, stability2$selprop)
      stability1$Beta <- abind::abind(stability1$Beta, stability2$Beta, along = 1)
    }
    if (stability1$methods$type == "graphical_model") {
      stability1$selprop <- abind::abind(stability1$selprop, stability2$selprop, along = 3)
    }
  }

  if (order_output) {
    ids <- sort.list(stability1$Q)
    stability1$S <- stability1$S[ids, , drop = FALSE]
    stability1$Lambda <- stability1$Lambda[ids, , drop = FALSE]
    stability1$Q <- stability1$Q[ids, , drop = FALSE]
    stability1$Q_s <- stability1$Q_s[ids, , drop = FALSE]
    stability1$P <- stability1$P[ids, , drop = FALSE]
    stability1$PFER <- stability1$PFER[ids, , drop = FALSE]
    stability1$FDP <- stability1$FDP[ids, , drop = FALSE]
    stability1$S_2d <- stability1$S_2d[ids, , drop = FALSE]
    stability1$PFER_2d <- stability1$PFER_2d[ids, , drop = FALSE]
    stability1$FDP_2d <- stability1$FDP_2d[ids, , drop = FALSE]
    if (stability1$methods$type == "variable_selection") {
      stability1$selprop <- stability1$selprop[ids, , drop = FALSE]
      stability1$Beta <- stability1$Beta[ids, , , drop = FALSE]
    }
    if (stability1$methods$type == "graphical_model") {
      stability1$selprop <- stability1$selprop[, , ids, drop = FALSE]
    }
  }

  return(stability1)
}

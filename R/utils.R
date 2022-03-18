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
  if ("matrix" %in% class(object)) {
    theta <- object
  } else {
    # From igraph object
    if (class(object) == "igraph") {
      theta <- as.matrix(igraph::get.adjacency(object))
    }

    # From output of VariableSelection()
    if (class(object) == "variable_selection") {
      theta <- cbind(SelectedVariables(object))
      if (ncol(theta) == 1) {
        colnames(theta) <- "outcome"
      }
      theta <- Square(theta)
    }

    # From output of GraphicalModel()
    if (class(object) == "graphical_model") {
      theta <- Adjacency(object)
    }

    # From output of BiSelection()
    if (class(object) == "bi_selection") {
      if ("selected" %in% names(object)) {
        selected <- object$selected # PLS
      } else {
        selected <- object$selectedX # PCA
      }
      theta <- Square(selected)
    }

    # From output of SimulateRegression() or SimulateComponents()
    if (class(object) %in% c("simulation_regression", "simulation_components")) {
      theta <- cbind(object$theta)
      if (ncol(theta) == 1) {
        colnames(theta) <- "outcome"
      }
      theta <- Square(theta)
    }

    # From output of SimulateGraphical()
    if (class(object) == "simulation_graphical_model") {
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
#' \dontrun{
#' # Simulated links between two sets
#' set.seed(1)
#' mat <- matrix(sample(c(0, 1), size = 15, replace = TRUE),
#'   nrow = 5, ncol = 3
#' )
#'
#' # Adjacency matrix of a bipartite graph
#' Square(mat)
#' }
#'
#' @export
Square <- function(x) {
  # Assigning row and column names
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("row", 1:nrow(x))
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("col", 1:ncol(x))
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
#' @examples
#' \dontrun{
#' # Simulated categorical variable
#' cat <- c(rep(1, 3), rep(2, 5), rep(3, 2))
#'
#' # Dummy variables
#' dummy <- as.matrix(model.matrix(~ as.factor(cat) - 1))
#'
#' # Back to categories
#' DummyToCategories(dummy)
#' }
#'
#' @export
DummyToCategories <- function(x, verbose = FALSE) {
  x_original <- x
  x <- matrix(0, nrow = nrow(x_original), ncol = ncol(x_original))
  for (j in 1:ncol(x)) {
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


#' Matching arguments
#'
#' Returns a vector of overlapping character strings between \code{extra_args}
#' and arguments from function \code{FUN}. If \code{FUN} is taking \code{...} as
#' input, this function returns \code{extra_args}.
#'
#' @param extra_args vector of character strings.
#' @param FUN function.
#'
#' @return A vector of overlapping arguments.
#'
#' @examples
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   MatchingArguments(
#'     extra_args = list(scale = TRUE, lambda = 1),
#'     FUN = sgPLS::sPLS
#'   )
#' }
#' @export
MatchingArguments <- function(extra_args, FUN) {
  if ("..." %in% names(formals(FUN))) {
    out <- extra_args
  } else {
    ids <- which(names(extra_args) %in% names(formals(FUN)))
    out <- extra_args[ids]
  }
  return(out)
}

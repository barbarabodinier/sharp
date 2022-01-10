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


#' Pairwise co-membership
#'
#' Generates a symmetric and binary matrix indicating if a pair of features
#' belongs to the same cluster.
#'
#' @param groups vector of group membership.
#'
#' @return A symmetric and binary matrix.
#'
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
  diag(comembership) <- 0
  rownames(comembership) <- colnames(comembership) <- names(groups)

  return(comembership)
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

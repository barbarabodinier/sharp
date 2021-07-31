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

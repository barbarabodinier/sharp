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
  comembership <- matrix(0, nrow = length(groups), ncol = length(groups))
  for (i in 1:(nrow(comembership) - 1)) {
    for (j in (i + 1):ncol(comembership)) {
      comembership[i, j] <- ifelse(groups[i] == groups[j], yes = 1, no = 0)
    }
  }
  comembership <- comembership + t(comembership)
  return(comembership)
}

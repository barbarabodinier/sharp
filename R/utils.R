#' Transforms NA into NULL
#'
#' Returns a vector with no missing values or NULL if there are no non-missing
#' values.
#'
#' @param x input vector.
#'
#' @return a vector without missing values or NULL.
#'
#' @export
NAToNULL <- function(x) {
  if (any(!is.na(x))) {
    return(x = x[!is.na(x)])
  } else {
    return(NULL)
  }
}

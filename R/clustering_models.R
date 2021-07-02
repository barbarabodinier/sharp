#' Hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. This function is not using stability.
#'
#' @inheritParams Clustering
#' @param ... additional parameters passed to \code{\link[stats]{hclust}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' myhclust <- HierarchicalClustering(xdata = t(simul$data), Lambda = 1:20)
#' @export
HierarchicalClustering <- function(xdata, Lambda = NULL, scale = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing for clustering of columns
  xdata <- t(xdata)

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Extracting relevant extra arguments (distance)
  ids <- which(names(extra_args) %in% names(formals(stats::dist)))
  ids <- ids[!ids %in% c("x")]

  # Computing pairwise distances
  mydistance <- do.call(stats::dist, args = c(list(x = xdata), extra_args[ids]))

  # Extracting relevant extra arguments (hclust)
  ids <- which(names(extra_args) %in% names(formals(stats::hclust)))
  ids <- ids[!ids %in% c("d")]

  # Running hierarchical clustering
  myclust <- do.call(stats::hclust, args = c(list(d = mydistance), extra_args[ids]))

  # Initialisation of array storing co-membership matrices
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(Lambda)))

  # Defining clusters
  mygroups <- do.call(stats::cutree, args = list(tree = myclust, k = Lambda))
  for (i in 1:nrow(Lambda)) {
    adjacency[, , i] <- CoMembership(groups = mygroups[, i])
  }

  return(adjacency)
}

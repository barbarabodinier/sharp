#' Hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. This function is not using stability.
#'
#' @inheritParams Clustering
#' @param rows logical indicating if clusters of rows (\code{TRUE}) or columns (\code{FALSE})
#'   should be conducted.
#' @param ... additional parameters passed to \code{\link[stats]{dist}} or
#'   \code{\link[stats]{hclust}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Hierarchical clustering
#' myhclust <- HierarchicalClustering(xdata = simul$data, Lambda = 1:20)
#' @export
HierarchicalClustering <- function(xdata, Lambda = NULL, scale = TRUE, rows = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

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
  if (is.null(dim(mygroups))) {
    mygroups <- cbind(mygroups)
  }
  for (i in 1:nrow(Lambda)) {
    adjacency[, , i] <- CoMembership(groups = mygroups[, i])
  }

  return(adjacency)
}


#' K-means clustering
#'
#' Runs k-means clustering using implementation from
#' \code{\link[stats]{kmeans}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[stats]{kmeans}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # k-means clustering
#' mykmeans <- KMeansClustering(xdata = simul$data, Lambda = 1:20)
#' @export
KMeansClustering <- function(xdata, Lambda = NULL, scale = TRUE, rows = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(Lambda)))

  # Extracting relevant extra arguments (kmeans)
  ids <- which(names(extra_args) %in% names(formals(stats::kmeans)))
  ids <- ids[!ids %in% c("x", "centers")]

  # Running k-means clustering
  for (k in 1:nrow(Lambda)) {
    if (Lambda[k, 1] < nrow(xdata)) {
      myclust <- do.call(stats::kmeans, args = c(list(x = xdata, centers = Lambda[k, 1]), extra_args[ids]))
      mygroups <- myclust$cluster
      adjacency[, , k] <- CoMembership(groups = mygroups)
    }
  }

  return(adjacency)
}


#' Partitioning Around Medoids
#'
#' Runs Partitioning Around Medoids (PAM) clustering using implementation from
#' \code{\link[cluster]{pam}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[stats]{dist}} or
#'   \code{\link[cluster]{pam}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # PAM clustering
#' mypam <- PAMClustering(xdata = simul$data, Lambda = 1:20)
#' @export
PAMClustering <- function(xdata, Lambda = NULL, scale = TRUE, rows = TRUE, ...) {
  # Checking cluster package is installed
  if (!requireNamespace("cluster")) {
    stop("This function requires the 'cluster' package.")
  }
  
  # Storing extra arguments
  extra_args <- list(...)

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(Lambda)))

  # Extracting relevant extra arguments (distance)
  ids <- which(names(extra_args) %in% names(formals(stats::dist)))
  ids <- ids[!ids %in% c("x")]

  # Computing pairwise distances
  mydistance <- do.call(stats::dist, args = c(list(x = xdata), extra_args[ids]))

  # Extracting relevant extra arguments (pam)
  ids <- which(names(extra_args) %in% names(formals(cluster::pam)))
  ids <- ids[!ids %in% c("x", "k", "diss", "cluster.only")]

  # Running k-means clustering
  for (k in 1:nrow(Lambda)) {
    if (Lambda[k, 1] < nrow(xdata)) {
      mygroups <- do.call(cluster::pam, args = c(
        list(x = mydistance, k = Lambda[k, 1], diss = TRUE, cluster.only = TRUE),
        extra_args[ids]
      ))
      adjacency[, , k] <- CoMembership(groups = mygroups)
    }
  }

  return(adjacency)
}

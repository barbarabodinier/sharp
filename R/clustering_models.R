#' Hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. This function is not using stability.
#'
#' @inheritParams Clustering
#' @param rows logical indicating if clusters of rows (\code{TRUE}) or columns (\code{FALSE})
#'   should be inferred.
#' @param ... additional parameters passed to \code{\link[stats]{dist}} or
#'   \code{\link[stats]{hclust}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @family clustering algorithms
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

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
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
#' @family clustering algorithms
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

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
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
    if (Lambda[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (Lambda[k, 1] < nrow(xdata)) {
        myclust <- do.call(stats::kmeans, args = c(list(x = xdata, centers = Lambda[k, 1]), extra_args[ids]))
        mygroups <- myclust$cluster
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(adjacency)
}


#' Model-based clustering
#'
#' Runs clustering with Gaussian Mixture Models (GMM) using implementation from
#' \code{\link[mclust]{Mclust}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[mclust]{Mclust}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @family clustering algorithms
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Clustering using Gaussian Mixture Models
#' mygmm <- GMMClustering(xdata = simul$data, Lambda = 1:20)
#' @importFrom mclust mclustBIC
#' @export
GMMClustering <- function(xdata, Lambda = NULL, scale = TRUE, rows = TRUE, ...) {
  # Checking mclust package is installed
  if (!requireNamespace("mclust")) {
    stop("This function requires the 'mclust' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(Lambda)))

  # Extracting relevant extra arguments (kmeans)
  ids <- which(names(extra_args) %in% names(formals(stats::kmeans)))
  ids <- ids[!ids %in% c("data", "G", "verbose")]

  # Running k-means clustering
  for (k in 1:nrow(Lambda)) {
    if (Lambda[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (Lambda[k, 1] < nrow(xdata)) {
        myclust <- do.call(mclust::Mclust, args = c(list(data = xdata, G = Lambda[k, 1], verbose = FALSE), extra_args[ids]))
        mygroups <- myclust$classification
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(adjacency)
}


#' Partitioning Around Medoids
#'
#' Runs Partitioning Around Medoids (PAM) clustering using implementation from
#' \code{\link[cluster]{pam}}. This is also known as the k-medoids algorithm.
#' This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[stats]{dist}} or
#'   \code{\link[cluster]{pam}}.
#'
#' @details Faster implementations of the algorithm can be chosen via the
#'   parameter \code{pamonce} (see \code{\link[cluster]{pam}}).
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @family clustering algorithms
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

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
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
    if (Lambda[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (Lambda[k, 1] < nrow(xdata)) {
        mygroups <- do.call(cluster::pam, args = c(
          list(x = mydistance, k = Lambda[k, 1], diss = TRUE, cluster.only = TRUE),
          extra_args[ids]
        ))
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(adjacency)
}


#' Clustering Large Applications
#'
#' Runs Clustering Large Applications (CLARA) algorithm using implementation
#' from \code{\link[cluster]{clara}}. This function applies the Partitioning
#' Around Medoids (PAM) algorithm on subsamples of the features to cluster. It
#' makes it more efficient for high-dimensional datasets. This function is not
#' using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[cluster]{clara}}.
#'
#' @return An array with binary and symmetric co-membership matrices.
#'
#' @family clustering algorithms
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # PAM clustering
#' mypam <- CLARAClustering(xdata = simul$data, Lambda = 1:20)
#' @export
CLARAClustering <- function(xdata, Lambda = NULL, scale = TRUE, rows = TRUE, ...) {
  # Checking cluster package is installed
  if (!requireNamespace("cluster")) {
    stop("This function requires the 'cluster' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Transposing for clustering of columns
  if (!rows) {
    xdata <- t(xdata)
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(Lambda)))

  # Extracting relevant extra arguments (clara)
  ids <- which(names(extra_args) %in% names(formals(cluster::clara)))
  ids <- ids[!ids %in% c("x", "k", "cluster.only", "stand")]

  # Running k-means clustering
  for (k in 1:nrow(Lambda)) {
    if (Lambda[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (Lambda[k, 1] < nrow(xdata)) {
        mygroups <- do.call(cluster::clara, args = c(
          list(x = xdata, k = Lambda[k, 1], stand = scale, cluster.only = TRUE),
          extra_args[ids]
        ))
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(adjacency)
}

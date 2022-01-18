#' Hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. This function is not using stability.
#'
#' @param xdata data matrix with observations as rows and variables as columns.
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc} is
#'   not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param rows logical indicating if clusters of rows (\code{TRUE}) or columns (\code{FALSE})
#'   should be inferred.
#' @param scale logical indicating if the data should be scaled to ensure that
#'   all variables contribute equally to the clustering of the observations.
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
#' myhclust <- HierarchicalClustering(xdata = simul$data, nc = 1:20)
#' @export
HierarchicalClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Extracting relevant extra arguments (dist)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x")]

  # Computing pairwise distances
  mydistance <- do.call(stats::dist, args = c(list(x = xdata), tmp_extra_args))

  # Extracting relevant extra arguments (hclust)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::hclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("d")]

  # Running hierarchical clustering
  myclust <- do.call(stats::hclust, args = c(list(d = mydistance), tmp_extra_args))

  # Initialisation of array storing co-membership matrices
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Defining clusters
  mygroups <- do.call(stats::cutree, args = list(tree = myclust, k = nc))
  if (is.null(dim(mygroups))) {
    mygroups <- cbind(mygroups)
  }
  for (i in 1:nrow(nc)) {
    adjacency[, , i] <- CoMembership(groups = mygroups[, i])
  }

  return(list(comembership = adjacency))
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
#' mykmeans <- KMeansClustering(xdata = simul$data, nc = 1:20)
#' @export
KMeansClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::kmeans)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "centers")]

  # Running k-means clustering
  for (k in 1:nrow(nc)) {
    if (nc[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (nc[k, 1] < nrow(xdata)) {
        myclust <- do.call(stats::kmeans, args = c(list(x = xdata, centers = nc[k, 1]), tmp_extra_args))
        mygroups <- myclust$cluster
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(list(comembership = adjacency))
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
#' mygmm <- GMMClustering(xdata = simul$data, nc = 1:30)
#' @export
GMMClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = mclust::Mclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("data", "G", "verbose")]

  # Running k-means clustering
  for (k in 1:nrow(nc)) {
    if (nc[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (nc[k, 1] < nrow(xdata)) {
        myclust <- do.call(mclust::Mclust, args = c(list(data = xdata, G = nc[k, 1], verbose = FALSE), tmp_extra_args))
        mygroups <- myclust$classification
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(list(comembership = adjacency))
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
#' mypam <- PAMClustering(xdata = simul$data, nc = 1:20)
#' mypam$comembership[, , 5]
#'
#' # Using another distance metric: argument passed to stats::dist()
#' mypam <- PAMClustering(xdata = simul$data, nc = 1:20, method = "manhattan")
#' mypam$comembership[, , 5]
#'
#' # Using specific medoids: arguments passed to cluster::pam()
#' mypam <- PAMClustering(
#'   xdata = simul$data, nc = 5,
#'   medoids = 1:5, do.swap = FALSE
#' )
#' plot(Graph(
#'   adjacency = mypam$comembership[, , 1],
#'   node_colour = c(rep("red", 5), rep("blue", 15)),
#'   satellites = TRUE
#' ))
#' @export
PAMClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Extracting relevant extra arguments (dist)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x")]

  # Computing pairwise distances
  mydistance <- do.call(stats::dist, args = c(list(x = xdata), tmp_extra_args))

  # Extracting relevant extra arguments (pam)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = cluster::pam)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "k", "diss", "cluster.only", "stand")] # stand is ignored since x is dissimilarity matrix

  # Running k-means clustering
  for (k in 1:nrow(nc)) {
    if (nc[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (nc[k, 1] < nrow(xdata)) {
        mygroups <- do.call(cluster::pam, args = c(
          list(x = mydistance, k = nc[k, 1], diss = TRUE, cluster.only = TRUE),
          tmp_extra_args
        ))
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(list(comembership = adjacency))
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
#' mypam <- CLARAClustering(xdata = simul$data, nc = 1:20)
#' @export
CLARAClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = cluster::clara)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "k", "stand", "cluster.only")]

  # Running k-means clustering
  for (k in 1:nrow(nc)) {
    if (nc[k, 1] == 1) {
      adjacency[, , k] <- CoMembership(groups = rep(1, nrow(xdata)))
    } else {
      if (nc[k, 1] < nrow(xdata)) {
        mygroups <- do.call(cluster::clara, args = c(
          list(x = xdata, k = nc[k, 1], stand = scale, cluster.only = TRUE),
          tmp_extra_args
        ))
        adjacency[, , k] <- CoMembership(groups = mygroups)
      }
    }
  }

  return(list(comembership = adjacency))
}


#' Density-Based Spatial Clustering of Applications with Noise
#'
#' Runs Density-Based Spatial Clustering of Applications with Noise (DBSCAN)
#' algorithm using implementation from \code{\link[dbscan]{dbscan}}. This
#' function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[dbscan]{dbscan}}.
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
#' mydbscan <- DBSCANClustering(
#'   xdata = simul$data,
#'   nc = seq(6, 10, by = 0.1)
#' )
#' @export
DBSCANClustering <- function(xdata, nc = NULL, scale = TRUE, rows = TRUE, ...) {
  # Checking dbscan package is installed
  if (!requireNamespace("dbscan")) {
    stop("This function requires the 'dbscan' package.")
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

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = dbscan::dbscan)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "eps", "minPts")]

  # Running k-means clustering
  for (k in 1:nrow(nc)) {
    mygroups <- do.call(dbscan::dbscan, args = c(
      list(x = xdata, eps = nc[k, 1], minPts = 1),
      tmp_extra_args
    ))$cluster
    adjacency[, , k] <- CoMembership(groups = mygroups)
  }

  return(list(comembership = adjacency))
}

#' (Weighted) hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. If \code{Lambda} is provided, clustering is
#' applied on the weighted distance matrix calculated using the
#' \code{\link[rCOSA]{cosa2}} algorithm. Otherwise, distances are calculated
#' using \code{\link[stats]{dist}}. This function is not using stability.
#'
#' @inheritParams Clustering
#' @param Lambda vector of penalty parameters (see argument \code{lambda} in
#'   \code{\link[rCOSA]{cosa2}}). Unweighted distance matrices are used if
#'   \code{Lambda=NULL}.
#' @param distance character string indicating the type of distance to use. If
#'   \code{Lambda=NULL}, possible values include \code{"euclidean"},
#'   \code{"maximum"}, \code{"canberra"}, \code{"binary"}, and
#'   \code{"minkowski"} (see argument \code{method} in
#'   \code{\link[stats]{dist}}).  Otherwise, possible values include
#'   \code{"euclidean"} (\code{pwr=2}) or \code{"absolute"} (\code{pwr=1}) (see
#'   argument \code{pwr} in \code{\link[rCOSA]{cosa2}}).
#' @param ... additional parameters passed to \code{\link[stats]{hclust}},
#'   \code{\link[stats]{dist}}, or \code{\link[rCOSA]{cosa2}}. Parameters
#'   \code{niter} (default to 1) and \code{noit} (default to 100) correspond to
#'   the number of iterations in \code{\link[rCOSA]{cosa2}} to calculate weights
#'   and may need to be modified. Argument \code{pwr} in
#'   \code{\link[rCOSA]{cosa2}} is ignored, please provide \code{distance}
#'   instead.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @references \insertRef{rCOSA}{sharp}
#'
#'   \insertRef{COSA}{sharp}
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Hierarchical clustering
#' myhclust <- HierarchicalClustering(
#'   xdata = simul$data,
#'   nc = seq_len(20)
#' )
#'
#' # Weighted Hierarchical clustering (using COSA)
#' if (requireNamespace("rCOSA", quietly = TRUE)) {
#'   myhclust <- HierarchicalClustering(
#'     xdata = simul$data,
#'     weighted = TRUE,
#'     nc = seq_len(20),
#'     Lambda = c(0.2, 0.5)
#'   )
#' }
#' @export
HierarchicalClustering <- function(xdata, nc = NULL, Lambda = NULL,
                                   distance = "euclidean",
                                   linkage = "complete",
                                   ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Checking that pwr is not provided (rCOSA)
  if ("pwr" %in% names(extra_args)) {
    warning("Argument 'pwr' is ignored. Please provide 'distance' instead: 'euclidean' is equivalent to 'pwr = 2' and 'absolute' is equivalent to 'pwr = 1'.")
  }

  # Defining lX argument for cosa2
  if ("lX" %in% names(extra_args)) {
    lX <- extra_args["lX"]
  } else {
    lX <- rep(1, ncol(xdata))
  }

  # Preparing binary indicator for weighted clustering
  if (is.null(Lambda)) {
    weighted <- FALSE
    Lambda <- 0
  } else {
    weighted <- TRUE

    # Checking rCOSA package is installed
    CheckPackageInstalled("rCOSA")
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidean", "absolute")) {
      distance <- ifelse(distance == "euclidean", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidean' (L2-norm) or 'absolute' (L1-norm). The 'euclidean' distance was used.")
      distance <- 2
    }
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  if (weighted) {
    # Extracting relevant extra arguments (cosa2)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
    tmp_extra_args_cosa2 <- tmp_extra_args
  } else {
    # Extracting relevant extra arguments (dist)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
    tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]
  }

  # Extracting relevant extra arguments (hclust)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::hclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("d", "method")]
  tmp_extra_args_hclust <- tmp_extra_args

  # Initialisation of arrays storing co-membership matrices and weights
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * nrow(Lambda)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Making xdata a data frame for cosa2
  xdata <- as.data.frame(xdata)

  # Iterating over the pair of parameters
  id <- 0
  for (i in seq_len(nrow(Lambda))) {
    if (weighted) {
      # Computing weighted pairwise distances
      utils::capture.output({
        out <- do.call(rCOSA::cosa2, args = c(
          list(
            X = xdata,
            lX = lX,
            lambda = Lambda[i, 1],
            stand = 0,
            pwr = distance
          ),
          tmp_extra_args_cosa2
        ))
      })
      mydistance <- out$D
    } else {
      # Computing pairwise distances
      mydistance <- do.call(stats::dist, args = c(
        list(x = xdata, method = distance),
        tmp_extra_args_dist
      ))
    }

    # Running hierarchical clustering
    myclust <- do.call(stats::hclust, args = c(
      list(
        d = mydistance,
        method = linkage
      ),
      tmp_extra_args_hclust
    ))

    # Defining clusters
    mygroups <- do.call(stats::cutree, args = list(tree = myclust, k = nc))
    if (is.null(dim(mygroups))) {
      mygroups <- cbind(mygroups)
    }
    for (j in seq_len(nrow(nc))) {
      adjacency[, , id + j] <- CoMembership(groups = mygroups[, j])
      # weight[,,id + j] <- out$W
      if (weighted) {
        weight[id + j, ] <- apply(out$W, 2, stats::median)
      }
    }
    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
}


#' (Weighted) Partitioning Around Medoids
#'
#' Runs Partitioning Around Medoids (PAM) clustering using implementation from
#' \code{\link[cluster]{pam}}. This is also known as the k-medoids algorithm. If
#' \code{Lambda} is provided, clustering is applied on the weighted distance
#' matrix calculated using the COSA algorithm as implemented in
#' \code{\link[rCOSA]{cosa2}}. Otherwise, distances are calculated using
#' \code{\link[stats]{dist}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[cluster]{pam}},
#'   \code{\link[stats]{dist}}, or \code{\link[rCOSA]{cosa2}}. If
#'   \code{weighted=TRUE}, parameters \code{niter} (default to 1) and
#'   \code{noit} (default to 100) correspond to the number of iterations in
#'   \code{\link[rCOSA]{cosa2}} to calculate weights and may need to be
#'   modified.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @references \insertRef{rCOSA}{sharp}
#'
#'   \insertRef{COSA}{sharp}
#'
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#'   # PAM clustering
#'   myclust <- PAMClustering(
#'     xdata = simul$data,
#'     nc = seq_len(20)
#'   )
#'
#'   # Weighted PAM clustering (using COSA)
#'   if (requireNamespace("rCOSA", quietly = TRUE)) {
#'     myclust <- PAMClustering(
#'       xdata = simul$data,
#'       nc = seq_len(20),
#'       Lambda = c(0.2, 0.5)
#'     )
#'   }
#' }
#' @export
PAMClustering <- function(xdata, nc = NULL, Lambda = NULL,
                          distance = "euclidean",
                          ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Defining lX argument for cosa2
  if ("lX" %in% names(extra_args)) {
    lX <- extra_args["lX"]
  } else {
    lX <- rep(1, ncol(xdata))
  }

  # Preparing binary indicator for weighted clustering
  if (is.null(Lambda)) {
    weighted <- FALSE
    Lambda <- 0
  } else {
    weighted <- TRUE

    # Checking rCOSA package is installed
    CheckPackageInstalled("rCOSA")
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidean", "absolute")) {
      distance <- ifelse(distance == "euclidean", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidean' (L2-norm) or 'absolute' (L1-norm). The 'euclidean' distance was used.")
      distance <- 2
    }
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  if (weighted) {
    # Extracting relevant extra arguments (cosa2)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
    tmp_extra_args_cosa2 <- tmp_extra_args
  } else {
    # Extracting relevant extra arguments (dist)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
    tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]
  }

  # Extracting relevant extra arguments (pam)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = cluster::pam)
  tmp_extra_args_pam <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "k", "diss", "cluster.only", "stand")] # stand is ignored since x is dissimilarity matrix

  # Initialisation of arrays storing co-membership matrices and weights
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * nrow(Lambda)))
  # weight <- array(NA, dim = c(nrow(xdata), ncol(xdata), nrow(nc) * nrow(Lambda)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Making xdata a data frame for cosa2
  xdata <- as.data.frame(xdata)

  # Iterating over the pair of parameters
  id <- 0
  for (i in seq_len(nrow(Lambda))) {
    if (weighted) {
      # Computing weighted pairwise distances
      utils::capture.output({
        out <- do.call(rCOSA::cosa2, args = c(
          list(
            X = xdata,
            lX = lX,
            lambda = Lambda[i, 1],
            stand = 0,
            pwr = distance
          ),
          tmp_extra_args_cosa2
        ))
      })
      mydistance <- out$D
    } else {
      # Computing pairwise distances
      mydistance <- do.call(stats::dist, args = c(
        list(x = xdata, method = distance),
        tmp_extra_args_dist
      ))
    }

    # Running PAM clustering
    for (k in seq_len(nrow(nc))) {
      if (nc[k, 1] == 1) {
        adjacency[, , id + k] <- CoMembership(groups = rep(1, nrow(xdata)))
      } else {
        if (nc[k, 1] < nrow(xdata)) {
          mygroups <- do.call(cluster::pam, args = c(
            list(x = mydistance, k = nc[k, 1], diss = TRUE, cluster.only = TRUE),
            tmp_extra_args_pam
          ))
          adjacency[, , id + k] <- CoMembership(groups = mygroups)
          if (weighted) {
            weight[id + k, ] <- apply(out$W, 2, stats::median)
          }
        }
      }
    }

    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
}


#' (Weighted) density-based clustering
#'
#' Runs Density-Based Spatial Clustering of Applications with Noise (DBSCAN)
#' clustering using implementation from \code{\link[dbscan]{dbscan}}. This is
#' also known as the k-medoids algorithm. If \code{Lambda} is provided,
#' clustering is applied on the weighted distance matrix calculated using the
#' COSA algorithm as implemented in \code{\link[rCOSA]{cosa2}}. Otherwise,
#' distances are calculated using \code{\link[stats]{dist}}. This function is
#' not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param eps radius in density-based clustering, see
#'   \code{\link[dbscan]{dbscan}}.
#' @param ... additional parameters passed to \code{\link[dbscan]{dbscan}}
#'   (except for \code{minPts} which is fixed to \code{2}),
#'   \code{\link[stats]{dist}}, or \code{\link[rCOSA]{cosa2}}. If
#'   \code{weighted=TRUE}, parameters \code{niter} (default to 1) and
#'   \code{noit} (default to 100) correspond to the number of iterations in
#'   \code{\link[rCOSA]{cosa2}} to calculate weights and may need to be
#'   modified.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @references \insertRef{rCOSA}{sharp}
#'
#'   \insertRef{COSA}{sharp}
#'
#' @examples
#' if (requireNamespace("dbscan", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'   plot(simul)
#'
#'   # DBSCAN clustering
#'   myclust <- DBSCANClustering(
#'     xdata = simul$data,
#'     eps = seq(0, 2 * sqrt(ncol(simul$data) - 1), by = 0.1)
#'   )
#'
#'   # Weighted PAM clustering (using COSA)
#'   if (requireNamespace("rCOSA", quietly = TRUE)) {
#'     myclust <- DBSCANClustering(
#'       xdata = simul$data,
#'       eps = c(0.25, 0.5, 0.75),
#'       Lambda = c(0.2, 0.5)
#'     )
#'   }
#' }
#' @export
DBSCANClustering <- function(xdata,
                             nc = NULL,
                             eps = NULL,
                             Lambda = NULL,
                             distance = "euclidean",
                             ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Defining lX argument for cosa2
  if ("lX" %in% names(extra_args)) {
    lX <- extra_args["lX"]
  } else {
    lX <- rep(1, ncol(xdata))
  }

  # Preparing binary indicator for weighted clustering
  if (is.null(Lambda)) {
    weighted <- FALSE
    Lambda <- 0
  } else {
    weighted <- TRUE

    # Checking rCOSA package is installed
    CheckPackageInstalled("rCOSA")
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidean", "absolute")) {
      distance <- ifelse(distance == "euclidean", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidean' (L2-norm) or 'absolute' (L1-norm). The 'euclidean' distance was used.")
      distance <- 2
    }
  }

  # Re-formatting eps
  eps <- cbind(eps)

  # Initialising nc
  nc <- cbind(rep(NA, length(eps)))

  if (weighted) {
    # Extracting relevant extra arguments (cosa2)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
    tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
    tmp_extra_args_cosa2 <- tmp_extra_args
  } else {
    # Extracting relevant extra arguments (dist)
    tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
    tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]
  }

  # Extracting relevant extra arguments (pam)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = dbscan::dbscan)
  tmp_extra_args_dbscan <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "eps", "minPts", "row")]

  # Initialisation of arrays storing co-membership matrices and weights
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * nrow(Lambda)))
  # weight <- array(NA, dim = c(nrow(xdata), ncol(xdata), nrow(nc) * nrow(Lambda)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))
  nc_full <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = 1)

  # Initialising noise status
  noise <- matrix(0, nrow = nrow(xdata), ncol = nrow(nc) * nrow(Lambda))
  rownames(noise) <- rownames(xdata)

  # Making xdata a data frame for cosa2
  xdata <- as.data.frame(xdata)

  # Iterating over the pair of parameters
  id <- 0
  for (i in seq_len(nrow(Lambda))) {
    if (weighted) {
      # Computing weighted pairwise distances
      utils::capture.output({
        out <- do.call(rCOSA::cosa2, args = c(
          list(
            X = xdata,
            lX = lX,
            lambda = Lambda[i, 1],
            stand = 0,
            pwr = distance
          ),
          tmp_extra_args_cosa2
        ))
      })
      mydistance <- out$D
    } else {
      # Computing pairwise distances
      mydistance <- do.call(stats::dist, args = c(
        list(x = xdata, method = distance),
        tmp_extra_args_dist
      ))
    }

    # Running DBSCAN clustering
    for (k in seq_len(nrow(nc))) {
      # Running algorithm
      mygroups <- do.call(dbscan::dbscan, args = c(
        list(x = stats::as.dist(mydistance), eps = eps[k, 1], minPts = 2),
        tmp_extra_args_dbscan
      ))$cluster
      nc_full[id + k, 1] <- length(unique(mygroups))

      # Filling in co-membership matrix
      if (nc_full[k, 1] == 1) {
        adjacency[, , id + k] <- CoMembership(groups = rep(1, nrow(xdata)))
      } else {
        adjacency[, , id + k] <- CoMembership(groups = mygroups + 1)
      }
      if (weighted) {
        weight[id + k, ] <- apply(out$W, 2, stats::median)
      }

      # Extracting noise variables
      if (any(mygroups == 0)) {
        noise[which(mygroups == 0), id + k] <- 1
      }
    }

    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight, noise = noise, nc = nc_full))
}


#' (Sparse) K-means clustering
#'
#' Runs k-means clustering using implementation from
#' \code{\link[stats]{kmeans}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param Lambda vector of penalty parameters (see argument \code{wbounds} in
#'   \code{\link[sparcl]{KMeansSparseCluster}}).
#' @param ... additional parameters passed to \code{\link[stats]{kmeans}} (if
#'   \code{Lambda} is \code{NULL}) or \code{\link[sparcl]{KMeansSparseCluster}}.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @references \insertRef{SparseClustering}{sharp}
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # K means clustering
#' mykmeans <- KMeansClustering(xdata = simul$data, nc = seq_len(20))
#'
#' # Sparse K means clustering
#' if (requireNamespace("sparcl", quietly = TRUE)) {
#'   mykmeans <- KMeansClustering(
#'     xdata = simul$data, nc = seq_len(20),
#'     Lambda = c(2, 5)
#'   )
#' }
#'
#' @export
KMeansClustering <- function(xdata, nc = NULL,
                             Lambda = NULL,
                             ...) {
  if (is.null(Lambda)) {
    return(UnweightedKMeansClustering(xdata = xdata, nc = nc, ...))
  } else {
    return(SparseKMeansClustering(xdata = xdata, nc = nc, Lambda = Lambda, ...))
  }
}


#' Model-based clustering
#'
#' Runs clustering with Gaussian Mixture Models (GMM) using implementation from
#' \code{\link[mclust]{Mclust}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[mclust]{Mclust}}.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Clustering using Gaussian Mixture Models
#' mygmm <- GMMClustering(xdata = simul$data, nc = seq_len(30))
#' @export
GMMClustering <- function(xdata, nc = NULL,
                          ...) {
  # Checking mclust package is installed
  if (!requireNamespace("mclust")) {
    stop("This function requires the 'mclust' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  Lambda <- cbind(0)
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = mclust::Mclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("data", "G", "verbose")]

  # Running gaussian mixture model
  for (k in seq_len(nrow(nc))) {
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

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
}


#' Unweighted K-means clustering
#'
#' Runs k-means clustering using implementation from
#' \code{\link[stats]{kmeans}}. This function is not using stability.
#'
#' @inheritParams HierarchicalClustering
#' @param ... additional parameters passed to \code{\link[stats]{kmeans}}.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @keywords internal
UnweightedKMeansClustering <- function(xdata, nc = NULL,
                                       ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Initialisation of array storing co-membership matrices
  Lambda <- cbind(0)
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::kmeans)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "centers")]

  # Running k-means clustering
  for (k in seq_len(nrow(nc))) {
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

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
}


#' Sparse K means clustering
#'
#' Runs sparse K means clustering using implementation from
#' \code{\link[sparcl]{KMeansSparseCluster}}. This function is not using
#' stability.
#'
#' @inheritParams HierarchicalClustering
#' @param Lambda vector of penalty parameters (see argument \code{wbounds} in
#'   \code{\link[sparcl]{KMeansSparseCluster}}).
#' @param ... additional parameters passed to
#'   \code{\link[sparcl]{KMeansSparseCluster}}.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @references \insertRef{SparseClustering}{sharp}
#'
#' @keywords internal
SparseKMeansClustering <- function(xdata, nc = NULL, Lambda, ...) {
  # Checking sparcl package is installed
  if (!requireNamespace("sparcl")) {
    stop("This function requires the 'sparcl' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    nc <- cbind(seq(1, nrow(xdata)))
  }

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sparcl::KMeansSparseCluster)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "wbound", "silent", "K")]

  # Defining the number of iterations
  n_iter_lambda <- nrow(Lambda)
  Lambda_iter <- Lambda

  # Initialisation of array storing co-membership matrices
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * n_iter_lambda))
  weight <- matrix(NA, nrow = nrow(nc) * n_iter_lambda, ncol = ncol(xdata))

  # Iterating over the pair of parameters
  id <- 0
  for (i in seq_len(n_iter_lambda)) {
    for (j in seq_len(nrow(nc))) {
      # Running sparse K means clustering
      myclust <- tryCatch(
        do.call(sparcl::KMeansSparseCluster, args = c(
          list(x = xdata, K = nc[j], wbound = Lambda[i, 1], silent = TRUE),
          tmp_extra_args
        )),
        error = function(e) {
          NULL
        }
      )

      # Defining clusters
      if (!is.null(myclust)) {
        mygroups <- myclust[[1]]$Cs
        adjacency[, , id + j] <- CoMembership(groups = mygroups)
        weight[id + j, ] <- myclust[[1]]$ws
      } else {
        # Single cluster and unique weights in case of error
        mygroups <- rep(1, nrow(xdata))
        adjacency[, , id + j] <- CoMembership(groups = mygroups)
        weight[id + j, ] <- rep(1, ncol(xdata))
      }
    }
    id <- id + nrow(nc)
  }

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(
    comembership = adjacency,
    weight = weight,
    Lambda = Lambda_iter
  ))
}

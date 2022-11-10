#' (Weighted) hierarchical clustering
#'
#' Runs hierarchical clustering using implementation from
#' \code{\link[stats]{hclust}}. If \code{Lambda} is provided, clustering is
#' applied on the weighted distance matrix calculated using the COSA algorithm
#' as implemented in \code{\link[rCOSA]{cosa2}}. Otherwise, distances are
#' calculated using \code{\link[stats]{dist}}. This function is not using
#' stability.
#'
#' @inheritParams Clustering
#' @param Lambda vector of penalty parameters (see argument \code{lambda} from
#'   \code{\link[rCOSA]{cosa2}}). Classical (unweighted) clustering is used if
#'   \code{Lambda=NULL}.
#' @param distance character string indicating the type of distance to use. If
#'   \code{weighted=FALSE}, possible values include \code{"euclidian"},
#'   \code{"maximum"}, \code{"canberra"}, \code{"binary"}, and
#'   \code{"minkowski"} (see argument \code{method} in
#'   \code{\link[stats]{dist}}).  If \code{weighted=TRUE}, possible values
#'   include \code{"euclidian"} (\code{pwr=2}) or \code{"absolute"}
#'   (\code{pwr=1}) (see argument \code{pwr} from \code{\link[rCOSA]{cosa2}}).
#' @param ... additional parameters passed to \code{\link[stats]{hclust}},
#'   \code{\link[stats]{dist}}, or \code{\link[rCOSA]{cosa2}}. Parameters
#'   \code{niter} (default to 1) and \code{noit} (default to 100) correspond to
#'   the number of iterations in \code{\link[rCOSA]{cosa2}} to calculate weights
#'   and may need to be modified.
#'
#' @return A list with: \item{comembership}{an array of binary and symmetric
#'   co-membership matrices.} \item{weights}{a matrix of median weights by
#'   feature.}
#'
#' @family clustering algorithms
#'
#' @references \insertRef{COSA}{sharp}
#'
#' \insertRef{rCOSA}{sharp}
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
#'   nc = 1:20
#' )
#'
#' # Weighted Hierarchical clustering (using COSA)
#' myhclust <- HierarchicalClustering(
#'   xdata = simul$data,
#'   weighted = TRUE,
#'   nc = 1:20,
#'   Lambda = c(0.2, 0.5)
#' )
#' @export
HierarchicalClustering <- function(xdata, nc = NULL, Lambda = NULL,
                                   distance = "euclidian", linkage = "complete",
                                   scale = TRUE, row = TRUE, ...) {
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
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidian", "absolute")) {
      distance <- ifelse(distance == "euclidian", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidian' (L2-norm) or 'absolute' (L1-norm). The 'euclidian' distance was used.")
      distance <- 2
    }
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Transposing if clustering of columns
  if (!row) {
    xdata <- t(xdata)
  }

  # Re-formatting nc
  if (!is.null(nc)) {
    if (is.vector(nc)) {
      nc <- cbind(nc)
    }
  } else {
    if (row) {
      nc <- cbind(seq(1, nrow(xdata)))
    } else {
      nc <- cbind(seq(1, ncol(xdata)))
    }
  }

  # Extracting relevant extra arguments (cosa2)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
  tmp_extra_args_cosa2 <- tmp_extra_args

  # Extracting relevant extra arguments (dist)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
  tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]

  # Extracting relevant extra arguments (hclust)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::hclust)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("d", "method")]
  tmp_extra_args_hclust <- tmp_extra_args

  # Initialisation of arrays storing co-membership matrices and weights
  adjacency <- array(NA, dim = c(nrow(xdata), nrow(xdata), nrow(nc) * nrow(Lambda)))
  # weight <- array(NA, dim = c(nrow(xdata), ncol(xdata), nrow(nc) * nrow(Lambda)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

  # Making xdata a data frame for cosa2
  xdata <- as.data.frame(xdata)

  # Iterating over the pair of parameters
  id <- 0
  for (i in 1:nrow(Lambda)) {
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
    for (j in 1:nrow(nc)) {
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
#' @references \insertRef{COSA}{sharp}
#'
#' \insertRef{rCOSA}{sharp}
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # PAM clustering
#' myclust <- PAMClustering(
#'   xdata = simul$data,
#'   nc = 1:20
#' )
#'
#' # Weighted PAM clustering (using COSA)
#' myclust <- PAMClustering(
#'   xdata = simul$data,
#'   nc = 1:20,
#'   Lambda = c(0.2, 0.5)
#' )
#' @export
PAMClustering <- function(xdata, nc = NULL, Lambda = NULL,
                          distance = "euclidian",
                          scale = TRUE, ...) {
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
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidian", "absolute")) {
      distance <- ifelse(distance == "euclidian", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidian' (L2-norm) or 'absolute' (L1-norm). The 'euclidian' distance was used.")
      distance <- 2
    }
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

  # Extracting relevant extra arguments (cosa2)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
  tmp_extra_args_cosa2 <- tmp_extra_args

  # Extracting relevant extra arguments (dist)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
  tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]

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
  for (i in 1:nrow(Lambda)) {
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
    for (k in 1:nrow(nc)) {
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
#' @references \insertRef{COSA}{sharp}
#'
#'   \insertRef{rCOSA}{sharp}
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#' plot(simul)
#'
#' # DBSCAN clustering
#' myclust <- DBSCANClustering(
#'   xdata = simul$data,
#'   eps = seq(0, 2 * sqrt(ncol(simul$data) - 1), by = 0.1)
#' )
#'
#' # Weighted PAM clustering (using COSA)
#' myclust <- DBSCANClustering(
#'   xdata = simul$data,
#'   eps = c(0.25, 0.5, 0.75),
#'   Lambda = c(0.2, 0.5)
#' )
#' @export
DBSCANClustering <- function(xdata,
                             nc = NULL,
                             eps = NULL,
                             Lambda = NULL,
                             distance = "euclidian",
                             scale = TRUE, ...) {
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
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    Lambda <- cbind(Lambda)
  }

  # Adapting distance to cosa2
  if (weighted) {
    if (distance %in% c("euclidian", "absolute")) {
      distance <- ifelse(distance == "euclidian", yes = 2, no = 1)
    } else {
      warning("Invalid input for argument 'distance'. For COSA clustering, possible values are: 'euclidian' (L2-norm) or 'absolute' (L1-norm). The 'euclidian' distance was used.")
      distance <- 2
    }
  }

  # Scaling the data
  if (scale) {
    xdata <- scale(xdata)
  }

  # Re-formatting eps
  eps <- cbind(eps)

  # Initialising nc
  nc <- cbind(rep(NA, length(eps)))

  # Extracting relevant extra arguments (cosa2)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rCOSA::cosa2)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "lambda", "stand", "pwr")]
  tmp_extra_args_cosa2 <- tmp_extra_args

  # Extracting relevant extra arguments (dist)
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::dist)
  tmp_extra_args_dist <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "method")]

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
  for (i in 1:nrow(Lambda)) {
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
    for (k in 1:nrow(nc)) {
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


#' K-means clustering
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
KMeansClustering <- function(xdata, nc = NULL, scale = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

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
  Lambda <- cbind(0)
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

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

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
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
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10), pk = 50)
#'
#' # Clustering using Gaussian Mixture Models
#' mygmm <- GMMClustering(xdata = simul$data, nc = 1:30)
#' @export
GMMClustering <- function(xdata, nc = NULL, scale = TRUE, ...) {
  # Checking mclust package is installed
  if (!requireNamespace("mclust")) {
    stop("This function requires the 'mclust' package.")
  }

  # Storing extra arguments
  extra_args <- list(...)

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
  Lambda <- cbind(0)
  adjacency <- array(0, dim = c(nrow(xdata), nrow(xdata), nrow(nc)))
  weight <- matrix(NA, nrow = nrow(nc) * nrow(Lambda), ncol = ncol(xdata))

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

  # Setting row and column names
  rownames(weight) <- paste0("s", seq(0, nrow(weight) - 1))
  colnames(weight) <- colnames(xdata)

  return(list(comembership = adjacency, weight = weight))
}

#' Consensus clustering
#'
#' Performs consensus (weighted) clustering. The underlying algorithm (e.g.
#' hierarchical clustering) is run with different number of clusters \code{nc}.
#' In consensus weighed clustering, weighted distances are calculated using the
#' \code{\link[rCOSA]{cosa2}} algorithm with different penalty parameters
#' \code{Lambda}. The hyper-parameters are calibrated by maximisation of the
#' consensus score.
#'
#' @inheritParams VariableSelection
#' @param xdata data matrix with observations as rows and variables as columns.
#' @param tau subsample size.
#' @param Lambda vector of penalty parameters for weighted distance calculation.
#'   Only used for distance-based clustering, including for example
#'   \code{implementation=HierarchicalClustering},
#'   \code{implementation=PAMClustering}, or
#'   \code{implementation=DBSCANClustering}.
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc} is
#'   not provided, it is set to \code{seq(1, tau*nrow(xdata))}.
#' @param eps radius in density-based clustering, see
#'   \code{\link[dbscan]{dbscan}}. Only used if
#'   \code{implementation=DBSCANClustering}.
#' @param implementation function to use for clustering. Possible functions
#'   include \code{\link{HierarchicalClustering}} (hierarchical clustering),
#'   \code{\link{PAMClustering}} (Partitioning Around Medoids),
#'   \code{\link{KMeansClustering}} (k-means) and \code{\link{GMMClustering}}
#'   (Gaussian Mixture Models). Alternatively, a user-defined function taking
#'   \code{xdata} and \code{Lambda} as arguments and returning a binary and
#'   symmetric matrix for which diagonal elements are equal to zero can be used.
#' @param scale logical indicating if the data should be scaled to ensure that
#'   all variables contribute equally to the clustering of the observations.
#' @param linkage character string indicating the type of linkage used in
#'   hierarchical clustering to define the stable clusters. Possible values
#'   include \code{"complete"}, \code{"single"} and \code{"average"} (see
#'   argument \code{"method"} in \code{\link[stats]{hclust}} for a full list).
#'   Only used if \code{implementation=HierarchicalClustering}.
#' @param row logical indicating if rows (if \code{row=TRUE}) or columns (if
#'   \code{row=FALSE}) contain the items to cluster.
#'
#' @details In consensus clustering, a clustering algorithm is applied on
#'   \code{K} subsamples of the observations with different numbers of clusters
#'   provided in \code{nc}. If \code{row=TRUE} (the default), the observations
#'   (rows) are the items to cluster. If \code{row=FALSE}, the variables
#'   (columns) are the items to cluster. For a given number of clusters, the
#'   consensus matrix \code{coprop} stores the proportion of iterations where
#'   two items were in the same estimated cluster, out of all iterations where
#'   both items were drawn in the subsample.
#'
#'   Stable cluster membership is obtained by applying a distance-based
#'   clustering method using \code{(1-coprop)} as distance (see
#'   \link{Clusters}).
#'
#'   These parameters can be calibrated by maximisation of a stability score
#'   (see \code{\link{ConsensusScore}}) calculated under the null hypothesis of
#'   equiprobability of co-membership.
#'
#'   It is strongly recommended to examine the calibration plot (see
#'   \code{\link{CalibrationPlot}}) to check that there is a clear maximum. The
#'   absence of a clear maximum suggests that the clustering is not stable,
#'   consensus clustering outputs should not be trusted in that case.
#'
#'   To ensure reproducibility of the results, the starting number of the random
#'   number generator is set to \code{seed}.
#'
#'   For parallelisation, stability selection with different sets of parameters
#'   can be run on \code{n_cores} cores. Using \code{n_cores > 1} creates a
#'   \code{\link[future]{multisession}}.
#'
#' @return An object of class \code{clustering}. A list with: \item{Sc}{a matrix
#'   of the best stability scores for different (sets of) parameters controlling
#'   the number of clusters and penalisation of attribute weights.} \item{nc}{a
#'   matrix of numbers of clusters.} \item{Lambda}{a matrix of regularisation
#'   parameters for attribute weights.} \item{Q}{a matrix of the average number
#'   of selected attributes by the underlying algorithm with different
#'   regularisation parameters.} \item{coprop}{an array of consensus matrices.
#'   Rows and columns correspond to items. Indices along the third dimension
#'   correspond to different parameters controlling the number of clusters and
#'   penalisation of attribute weights.} \item{selprop}{an array of selection
#'   proportions. Columns correspond to attributes. Rows correspond to different
#'   parameters controlling the number of clusters and penalisation of attribute
#'   weights.} \item{method}{a list with \code{type="clustering"} and values
#'   used for arguments \code{implementation}, \code{linkage}, and
#'   \code{resampling}.} \item{params}{a list with values used for arguments
#'   \code{K}, \code{tau}, \code{pk}, \code{n} (number of observations in
#'   \code{xdata}), and \code{seed}.} The rows of \code{Sc}, \code{nc},
#'   \code{Lambda}, \code{Q}, \code{selprop} and indices along the third
#'   dimension of \code{coprop} are ordered in the same way and correspond to
#'   parameter values stored in \code{nc} and \code{Lambda}.
#'
#' @family stability functions
#'
#' @seealso \code{\link{Resample}}, \code{\link{ConsensusScore}},
#'   \code{\link{HierarchicalClustering}}, \code{\link{PAMClustering}},
#'   \code{\link{KMeansClustering}}, \code{\link{GMMClustering}}
#'
#' @references \insertRef{OurConsensusClustering}{sharp}
#'
#'   \insertRef{rCOSA}{sharp}
#'
#'   \insertRef{COSA}{sharp}
#'
#'   \insertRef{ConsensusClustering}{sharp}
#'
#' @examples
#' \donttest{
#' # Consensus clustering
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(30, 30, 30), nu_xc = 1, ev_xc = 0.5
#' )
#' stab <- Clustering(xdata = simul$data)
#' print(stab)
#' CalibrationPlot(stab)
#' summary(stab)
#' Clusters(stab)
#' plot(stab)
#'
#' # Consensus weighted clustering
#' if (requireNamespace("rCOSA", quietly = TRUE)) {
#'   set.seed(1)
#'   simul <- SimulateClustering(
#'     n = c(30, 30, 30), pk = 20,
#'     theta_xc = c(rep(1, 10), rep(0, 10)),
#'     ev_xc = 0.9
#'   )
#'   stab <- Clustering(
#'     xdata = simul$data,
#'     Lambda = LambdaSequence(lmin = 0.1, lmax = 10, cardinal = 10),
#'     noit = 20, niter = 10
#'   )
#'   print(stab)
#'   CalibrationPlot(stab)
#'   summary(stab)
#'   Clusters(stab)
#'   plot(stab)
#'   WeightBoxplot(stab)
#' }
#' }
#' @export
Clustering <- function(xdata, nc = NULL, eps = NULL, Lambda = NULL,
                       K = 100, tau = 0.5, seed = 1, n_cat = 3,
                       implementation = HierarchicalClustering,
                       scale = TRUE,
                       linkage = "complete",
                       row = TRUE,
                       n_cores = 1, output_data = FALSE, verbose = TRUE, beep = NULL, ...) {
  # Visiting all possible numbers of clusters
  if (is.null(eps)) {
    if (is.null(nc)) {
      if (row) {
        nc <- cbind(seq(1, nrow(xdata) * tau * 0.5))
      } else {
        nc <- cbind(seq(1, ncol(xdata) * 0.5))
      }
    } else {
      if (row) {
        if (any(nc > (nrow(xdata) * tau))) {
          nc <- nc[nc <= (nrow(xdata) * tau)]
          if (length(nc) > 0) {
            warning(paste0(
              "Invalid input for argument 'nc'. The number of clusters can not exceed the subsample size: ",
              nrow(xdata) * tau, ". Invalid values have been removed."
            ))
            nc <- cbind(nc)
          } else {
            stop(paste0(
              "Invalid input for argument 'nc'. The number of clusters can not exceed the subsample size: ",
              nrow(xdata) * tau, "."
            ))
          }
        }
      }
    }
  } else {
    nc <- cbind(rep(NA, length(eps)))
  }

  # Error and warning messages
  resampling <- "subsampling" # only subsampling is available as bootstrap would give distance of zero
  PFER_method <- "MB"
  PFER_thr <- Inf # different interpretation in clustering
  FDP_thr <- Inf
  pk <- NULL
  lambda_other_blocks <- 0.1
  start <- "warm"
  Lambda_cardinal <- 50
  lambda_max <- NULL
  lambda_path_factor <- 0.001
  max_density <- 0.5
  pi_list <- seq(0.6, 0.9, by = 0.01)
  bigblocks <- bigblocks_vect <- blocks <- N_blocks <- nblocks <- PFER_thr_blocks <- FDP_thr_blocks <- NULL
  CheckInputClustering(
    xdata = xdata, Lambda = Lambda,
    pi_list = pi_list, K = K, tau = tau, seed = seed,
    implementation = implementation, scale = scale,
    resampling = resampling,
    verbose = verbose
  )

  # Stability selection and score
  if (n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
    mypar <- future.apply::future_lapply(X = seq_len(n_cores), future.seed = TRUE, FUN = function(k) {
      return(SerialClustering(
        xdata = xdata, Lambda = cbind(Lambda), nc = cbind(nc), eps = cbind(eps),
        K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
        implementation = implementation, scale = scale, linkage = linkage, row = row,
        output_data = output_data, verbose = FALSE, ...
      ))
    }) # keep pk for correct number of blocks etc
    future::plan(future::sequential)

    # Combining the outputs from parallel iterations
    out <- mypar[[1]]
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[i]], include_beta = TRUE))
    }
  } else {
    out <- SerialClustering(
      xdata = xdata, Lambda = cbind(Lambda), nc = cbind(nc), eps = cbind(eps),
      K = K, tau = tau, seed = seed, n_cat = n_cat,
      implementation = implementation, scale = scale, linkage = linkage, row = row,
      output_data = output_data, verbose = verbose, ...
    )
  }

  # Re-set the function names
  if ("methods" %in% names(out)) {
    myimplementation <- as.character(substitute(implementation))
    if (is.function(resampling)) {
      myresampling <- as.character(substitute(resampling))
    } else {
      myresampling <- resampling
    }
    out$methods$implementation <- myimplementation
    out$methods$resampling <- myresampling
    out$methods$type <- "clustering"
  }

  # Defining the class
  class(out) <- "clustering"

  # Making beep
  if (!is.null(beep)) {
    beepr::beep(sound = beep)
  }

  return(out)
}


#' Consensus clustering (internal)
#'
#' Performs consensus (weighted) clustering. The underlying algorithm (e.g.
#' hierarchical clustering) is run with different number of clusters \code{nc}.
#' In consensus weighed clustering, weighted distances are calculated using the
#' \code{\link[rCOSA]{cosa2}} algorithm with different penalty parameters
#' \code{Lambda}. The hyper-parameters are calibrated by maximisation of the
#' consensus score. This function uses a serial implementation and requires the grids of
#' hyper-parameters as input (for internal use only).
#'
#' @inheritParams Clustering
#'
#' @return A list with: \item{Sc}{a matrix
#'   of the best stability scores for different (sets of) parameters controlling
#'   the number of clusters and penalisation of attribute weights.} \item{nc}{a
#'   matrix of numbers of clusters.} \item{Lambda}{a matrix of regularisation
#'   parameters for attribute weights.} \item{Q}{a matrix of the average number
#'   of selected attributes by the underlying algorithm with different
#'   regularisation parameters.} \item{coprop}{an array of consensus matrices.
#'   Rows and columns correspond to items. Indices along the third dimension
#'   correspond to different parameters controlling the number of clusters and
#'   penalisation of attribute weights.} \item{selprop}{an array of selection
#'   proportions. Columns correspond to attributes. Rows correspond to different
#'   parameters controlling the number of clusters and penalisation of attribute
#'   weights.} \item{method}{a list with \code{type="clustering"} and values
#'   used for arguments \code{implementation}, \code{linkage}, and
#'   \code{resampling}.} \item{params}{a list with values used for arguments
#'   \code{K}, \code{tau}, \code{pk}, \code{n} (number of observations in
#'   \code{xdata}), and \code{seed}.} The rows of \code{Sc}, \code{nc},
#'   \code{Lambda}, \code{Q}, \code{selprop} and indices along the third
#'   dimension of \code{coprop} are ordered in the same way and correspond to
#'   parameter values stored in \code{nc} and \code{Lambda}.
#'
#' @keywords internal
SerialClustering <- function(xdata, nc, eps, Lambda,
                             K = 100, tau = 0.5, seed = 1, n_cat = 3,
                             implementation = HierarchicalClustering,
                             scale = TRUE, linkage = "complete", row = TRUE,
                             output_data = FALSE, verbose = TRUE, ...) {
  # Defining the full vectors of nc and Lambda
  if (is.null(Lambda)) {
    # nc_full <- cbind(nc)
    nc_full <- cbind(rep(0, length(nc)))
    Lambda_full <- cbind(rep(NA, nrow(nc)))
  } else {
    # nc_full <- cbind(rep(nc, nrow(Lambda)))
    nc_full <- cbind(rep(0, length(nc) * nrow(Lambda)))
    Lambda_full <- cbind(rep(Lambda[, 1], each = nrow(nc)))
  }
  rownames(Lambda_full) <- NULL

  # Defining resampling method (only subsampling is available as bootstrap would give distance of zero)
  resampling <- "subsampling"
  PFER_method <- "MB"
  PFER_thr <- Inf
  FDP_thr <- Inf

  # Initialising objects to be filled
  N <- N_block <- ncol(xdata)

  # Defining number of items
  if (row) {
    n_items <- nrow(xdata)
    item_names <- rownames(xdata)
    attribute_names <- colnames(xdata)
  } else {
    n_items <- ncol(xdata)
    item_names <- colnames(xdata)
    attribute_names <- rownames(xdata)
  }

  # Initialising array of co-membership proportions
  if (!is.null(Lambda)) {
    bigstab_obs <- array(0,
      dim = c(n_items, n_items, nrow(Lambda) * nrow(nc)),
      dimnames = list(item_names, item_names, NULL)
    )
  } else {
    bigstab_obs <- array(0,
      dim = c(n_items, n_items, nrow(nc)),
      dimnames = list(item_names, item_names, NULL)
    )
  }

  # Initialising for subsampling of the rows
  sampled_pairs <- matrix(0, nrow = nrow(xdata), ncol = nrow(xdata))
  rownames(sampled_pairs) <- colnames(sampled_pairs) <- rownames(xdata)

  # Initialising the array for contributing variables
  if (!is.null(Lambda)) {
    Beta <- array(0, dim = c(nrow(Lambda) * nrow(nc), ncol(xdata), K))
    rownames(Beta) <- paste0("s", seq(0, nrow(Beta) - 1))
    colnames(Beta) <- attribute_names
  } else {
    Beta <- NULL
  }

  # Initialising the proportion of noise variables
  bignoise <- matrix(0, nrow = n_items, ncol = nrow(nc_full))
  rownames(bignoise) <- item_names

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Computation of the selection proportions over Lambda
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  for (k in seq_len(K)) {
    # Subsampling observations
    s <- Resample(
      data = xdata, family = NULL, tau = tau, resampling = resampling, ...
    )
    Xsub <- xdata[s, , drop = FALSE]

    # Applying clustering algorithm
    mybeta <- ClusteringAlgo(
      xdata = Xsub,
      Lambda = Lambda, eps = eps, nc = nc, scale = scale, row = row,
      implementation = implementation, ...
    )

    # Summing the number of clusters
    nc_full <- nc_full + mybeta$nc

    # Summing noise variable status to get the number of iterations as noise
    bignoise[rownames(mybeta$noise), ] <- bignoise[rownames(mybeta$noise), ] + mybeta$noise

    # Storing weights, used to define set of selected variables
    if (!is.null(Lambda)) {
      Beta[rownames(mybeta$weight), colnames(mybeta$weight), k] <- mybeta$weight
    }

    # Storing co-membership status
    for (i in seq_len(dim(mybeta$comembership)[3])) {
      if (row) {
        bigstab_obs[s, s, i] <- bigstab_obs[s, s, i] + mybeta$comembership[, , i]
      } else {
        bigstab_obs[, , i] <- bigstab_obs[, , i] + mybeta$comembership[, , i]
      }
    }

    # Storing sampled pairs
    sampled_pairs[s, s] <- sampled_pairs[s, s] + 1

    if (verbose) {
      utils::setTxtProgressBar(pb, k / K)
    }
  }

  # Calculating the mean number of clusters over subsampling iterations
  nc_full <- nc_full / K

  # Computing the co-membership proportions
  for (i in seq_len(dim(bigstab_obs)[3])) {
    if (row) {
      bigstab_obs[, , i] <- bigstab_obs[, , i] / sampled_pairs
    } else {
      bigstab_obs[, , i] <- bigstab_obs[, , i] / K
    }
  }
  bigstab_obs[is.nan(bigstab_obs)] <- NA

  # Computing the noise proportion
  for (i in seq_len(nrow(bignoise))) {
    if (row) {
      bignoise[i, ] <- bignoise[i, ] / diag(sampled_pairs)[i]
    } else {
      bignoise[i, ] <- bignoise[i, ] / K
    }
  }

  if (verbose) {
    cat("\n")
  }

  # Imputation of missing values (accounting for the fact that 2 items potentially never get picked together)
  if (any(is.na(bigstab_obs))) {
    warning("Missing values in consensus matrix. These have been set to zero by default. Consider increasing the number of subsamples 'K'.")
    for (i in seq_len(dim(bigstab_obs)[3])) {
      if (any(is.na(bigstab_obs[, , i]))) {
        tmpmat <- bigstab_obs[, , i]
        tmpmat[which(is.na(tmpmat))] <- 0
        bigstab_obs[, , i] <- tmpmat
      }
    }
  }

  # Calibration of consensus clustering
  metrics2 <- matrix(NA, ncol = 1, nrow = dim(bigstab_obs)[3])
  for (i in seq_len(dim(bigstab_obs)[3])) {
    # Clustering on the consensus matrix
    sh_clust <- stats::hclust(stats::as.dist(1 - bigstab_obs[, , i]), method = linkage)

    # Identifying stable clusters
    theta <- CoMembership(groups = stats::cutree(sh_clust, k = ceiling(nc_full[i])))

    # Calculating the consensus score
    if (row) {
      metrics2[i] <- ConsensusScore(
        prop = (bigstab_obs[, , i])[upper.tri(bigstab_obs[, , i])],
        K = sampled_pairs[upper.tri(sampled_pairs)],
        theta = theta[upper.tri(theta)]
      )
    } else {
      metrics2[i] <- ConsensusScore(
        prop = (bigstab_obs[, , i])[upper.tri(bigstab_obs[, , i])],
        K = K,
        theta = theta[upper.tri(theta)]
      )
    }
  }

  # Computing the selection proportions
  if (!is.null(Lambda)) {
    bigstab_var <- matrix(NA, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab_var) <- colnames(Beta)
    rownames(bigstab_var) <- rownames(Beta)
    for (i in seq_len(nrow(Beta))) {
      for (j in seq_len(ncol(Beta))) {
        bigstab_var[i, j] <- sum(Beta[i, j, ] != 0) / K
      }
    }
    Q <- cbind(round(apply(bigstab_var, 1, sum)))
    rownames(Q) <- NULL
  } else {
    bigstab_var <- NULL
    Q <- cbind(rep(ncol(xdata), length(nc_full)))
  }

  if (verbose) {
    utils::setTxtProgressBar(pb, 1)
    cat("\n")
  }

  # Preparing outputs
  myimplementation <- as.character(substitute(implementation, env = parent.frame(n = 2)))
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  out <- list(
    Sc = metrics2,
    nc = nc_full,
    Lambda = Lambda_full,
    Q = Q,
    coprop = bigstab_obs,
    Beta = Beta,
    selprop = bigstab_var,
    sampled_pairs = sampled_pairs,
    bignoise = bignoise,
    # Q_s = metrics1$Q_s, P = metrics1$P,
    # PFER = metrics1$PFER, FDP = metrics1$FDP,
    # S_2d = metrics1$S_2d, PFER_2d = metrics1$PFER_2d, FDP_2d = metrics1$FDP_2d,
    methods = list(
      type = "clustering",
      implementation = myimplementation,
      linkage = linkage,
      resampling = myresampling
    ),
    params = list(
      K = K,
      tau = tau,
      pk = ncol(xdata),
      n = nrow(xdata),
      seed = seed
    )
  )

  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata))
  }

  # Defining the class
  class(out) <- "clustering"

  return(out)
}

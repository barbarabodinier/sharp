#' Consensus clustering
#'
#' Runs (sparse) consensus clustering models with different combinations of
#' parameters controlling the number of clusters in the underlying algorithm
#' (\code{nc}), penalty for the selection of variables contributing to the
#' clustering (\code{Lambda}) and thresholds co-membership proportions. These
#' three parameters are jointly calibrated by maximising the stability score of
#' the model. The threshold in selection proportions to identify the variables
#' driving the clustering is also chosen by maximising the stability score in
#' the calibrated clustering model. This function can be used to identify stable
#' clusters of observations sharing similar profiles, as well as variables
#' driving the clustering if a sparse clustering algorithm is used.
#'
#' @inheritParams VariableSelection
#' @param xdata data matrix with observations as rows and variables as columns.
#' @param Lambda vector of penalty parameters.
#' @param nc matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{nc} is
#'   not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param implementation function to use for clustering. Possible functions
#'   include \code{\link{HierarchicalClustering}} (hierarchical clustering),
#'   \code{\link{SparseHierarchicalClustering}} (sparse hierarchical
#'   clustering), \code{\link{KMeansClustering}} (k-means),
#'   \code{\link{GMMClustering}} (Gaussian Mixture Models),
#'   \code{\link{PAMClustering}} (Partioning Around Medoids),
#'   \code{\link{CLARAClustering}} (Clustering Large Applications) and
#'   \code{\link{DBSCANClustering}} (Density-Based Spatial Clustering of
#'   Applications with Noise). Alternatively, a user-defined function taking
#'   \code{xdata} and \code{Lambda} as arguments and returning a binary and
#'   symmetric matrix for which diagonal elements are equal to zero can be used.
#' @param scale logical indicating if the data should be scaled to ensure that
#'   all variables contribute equally to the clustering of the observations.
#'
#' @details To ensure reproducibility of the results, the state of the random
#'   number generator is fixed to \code{seed}. For parallelisation of the code,
#'   consensus clustering results produced with different \code{seed}s and all
#'   other parameters equal can be combined (more details in
#'   \code{\link{Combine}}).
#'
#' @return A list with: \item{S}{a matrix of the best (block-specific) stability
#'   scores for different parameters controlling the number of clusters in the
#'   underlying algorithm. } \item{Lambda}{a matrix of parameters controlling
#'   the number of clusters. } \item{Q}{a matrix of average numbers of
#'   co-members for different parameters controlling the number of clusters.}
#'   \item{Q_s}{a matrix of calibrated numbers of stable co-members for
#'   different parameters controlling the number of clusters in the underlying
#'   algorithm. } \item{P}{a matrix of calibrated thresholds in co-membership
#'   proportions for different parameters controlling the number of clusters in
#'   the underlying algorithm. } \item{PFER}{a matrix of the upper-bounds in
#'   PFER of calibrated consensus clustering models with different (sets of)
#'   parameters controlling the number of clusters in the underlying algorithm.}
#'   \item{FDP}{a matrix of the upper-bounds in FDP of calibrated consensus
#'   clustering models with different parameters controlling the number of
#'   clusters in the underlying algorithm.} \item{S_2d}{an array of stability
#'   scores obtained with different combinations of parameters. Columns
#'   correspond to different thresholds in co-membership proportions. }
#'   \item{PFER_2d}{an array of computed upper-bounds of PFER obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in co-membership proportions. } \item{FDP_2d}{an array of
#'   computed upper-bounds of FDP obtained with different combinations of
#'   parameters. Columns correspond to different thresholds in co-membership
#'   proportions. } \item{selprop}{an array of co-membership proportions. Rows
#'   and columns correspond to features being clustered (rows of \code{xdata}).
#'   Indices along the third dimension correspond to different parameters
#'   controlling the number of clusters in the underlying algorithm.}
#'   \item{methods}{a list with \code{type="clustering"}, \code{implementation},
#'   \code{resampling} and \code{PFER_method} values used for the run.}
#'   \item{param}{a list with values of other objects used for the run.} For all
#'   objects except \code{selprop} and those stored in \code{methods} or
#'   \code{params}, rows correspond to parameter values stored in the output
#'   \code{Lambda}.
#'
#' @family stability selection functions
#' @seealso \code{\link{Resample}}, \code{\link{StabilityScore}},
#'   \code{\link{SparseHierarchicalClustering}}, \code{\link{KMeansClustering}},
#'   \code{\link{GMMClustering}}, \code{\link{PAMClustering}}
#'
#' @references \insertRef{ourstabilityselection}{focus}
#'
#'   \insertRef{stabilityselectionMB}{focus}
#'
#'   \insertRef{stabilityselectionSS}{focus}
#'
#' @examples
#' \dontshow{
#' set.seed(1)
#' simul <- SimulateClustering(n = c(5, 5, 5), pk = 10)
#' stab <- Clustering(
#'   xdata = simul$data,
#'   nc = c(2, 3),
#'   Lambda = c(1.1, 1.2),
#'   K = 5,
#'   verbose = FALSE
#' )
#' mymembership <- Clusters(stab)
#' }
#' \dontrun{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(50, 50, 50), pk = 15,
#'   theta_xc = c(1, 1, 1, 1, 1, rep(0, 10)),
#'   ev = c(0.9, 0.8, 0.7, 0.7, 0.7, rep(0, 10)),
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = as.matrix(dist(simul$data)),
#'   colours = c("navy", "white", "red")
#' )
#' r2 <- NULL
#' for (k in 1:ncol(simul$data)) {
#'   mymodel <- lm(simul$data[, k] ~ as.factor(simul$theta))
#'   r2 <- c(r2, summary(mymodel)$r.squared)
#' }
#'
#' # Sparse consensus clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   Lambda = cbind(seq(1.1, 3, by = 0.05))
#' )
#'
#' # Stable clusters
#' CalibrationPlot(stab)
#' plot(Graph(Adjacency(stab), node_colour = simul$theta, satellites = TRUE))
#' ClusteringPerformance(theta = Clusters(stab), theta_star = simul$theta)
#'
#' # Stably contributing variables
#' SelectionProportions(stab)
#' SelectedVariables(stab)
#'
#' # Consensus clustering based on hierarchical clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = HierarchicalClustering
#' )
#' CalibrationPlot(stab, xlab = expression(italic(k)))
#' plot(Graph(Adjacency(stab), node_colour = simul$theta, satellites = TRUE))
#' ClusteringPerformance(theta = Clusters(stab), theta_star = simul$theta)
#'
#' # Consensus clustering based on k-means clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = KMeansClustering
#' )
#' CalibrationPlot(stab, xlab = expression(italic(k)))
#' plot(Graph(Adjacency(stab), node_colour = simul$theta, satellites = TRUE))
#' ClusteringPerformance(theta = Clusters(stab), theta_star = simul$theta)
#'
#' # Consensus clustering based on PAM clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = PAMClustering
#' )
#' CalibrationPlot(stab, xlab = expression(italic(k)))
#' plot(Graph(Adjacency(stab), node_colour = simul$theta, satellites = TRUE))
#' ClusteringPerformance(theta = Clusters(stab), theta_star = simul$theta)
#' }
#' @export
Clustering <- function(xdata, Lambda = NULL, nc = NULL,
                       pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3,
                       implementation = SparseHierarchicalClustering, scale = TRUE,
                       n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  # Visiting all possible numbers of clusters
  if (is.null(nc)) {
    nc <- cbind(seq(1, nrow(xdata) * tau))
  } else {
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

  # Setting fixed Lambda if not provided
  if (is.null(Lambda)) {
    Lambda <- 1.1
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
  bigblocks <- bigblocks_vect <- blocks <- N_blocks <- nblocks <- PFER_thr_blocks <- FDP_thr_blocks <- NULL
  CheckInputGraphical(
    xdata = xdata, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
    pi_list = pi_list, K = K, tau = tau, seed = seed, n_cat = n_cat,
    implementation = implementation, start = start, scale = scale,
    resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    Lambda_cardinal = Lambda_cardinal,
    lambda_max = lambda_max, lambda_path_factor = lambda_path_factor, max_density = max_density,
    verbose = verbose
  )

  # Check if parallelisation is possible (forking)
  if (.Platform$OS.type != "unix") {
    if (n_cores > 1) {
      warning("Invalid input for argument 'n_cores'. Parallelisation relies on forking, it is only available on Unix systems.")
    }
    n_cores <- 1
  }

  # Stability selection and score
  mypar <- parallel::mclapply(X = 1:n_cores, FUN = function(k) {
    return(SerialClustering(
      xdata = xdata, Lambda = cbind(Lambda), nc = cbind(nc),
      pi_list = pi_list, K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
      implementation = implementation,
      output_data = output_data, verbose = verbose, ...
    ))
  }) # keep pk for correct number of blocks etc

  # Combining the outputs from parallel iterations
  out <- mypar[[1]]
  if (n_cores > 1) {
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[2]], graph = TRUE))
    }
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

  return(out)
}


#' Consensus clustering (internal)
#'
#' Runs stability selection regression models with different combinations of
#' parameters controlling the sparsity of the underlying selection algorithm
#' (e.g. penalty parameter for regularised models) and thresholds in selection
#' proportions. These two parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features). This function uses a serial
#' implementation and requires the grid of parameters controlling the underlying
#' algorithm as input (for internal use only).
#'
#' @inheritParams Clustering
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   With \code{implementation="glmnet"}, \code{Lambda} contains penalty
#'   parameters.
#'
#' @return A list with: \item{S}{a matrix of the best stability scores for
#'   different parameters controlling the level of sparsity in the underlying
#'   algorithm.} \item{Lambda}{a matrix of parameters controlling the level of
#'   sparsity in the underlying algorithm.} \item{Q}{a matrix of the average
#'   number of selected features by underlying algorithm with different
#'   parameters controlling the level of sparsity.} \item{Q_s}{a matrix of the
#'   calibrated number of stably selected features with different parameters
#'   controlling the level of sparsity.} \item{P}{a matrix of calibrated
#'   thresholds in selection proportions for different parameters controlling
#'   the level of sparsity in the underlying algorithm.} \item{PFER}{a matrix of
#'   the upper-bounds in PFER of calibrated stability selection models with
#'   different parameters controlling the level of sparsity.} \item{FDP}{a
#'   matrix of the upper-bounds in FDP of calibrated stability selection models
#'   with different parameters controlling the level of sparsity.} \item{S_2d}{a
#'   matrix of stability scores obtained with different combinations of
#'   parameters. Columns correspond to different tresholds in selection
#'   proportions.} \item{selprop}{a matrix of selection proportions. Columns
#'   correspond to predictors from \code{xdata}.} \item{Beta}{an array of model
#'   coefficients. Columns correspond to predictors from \code{xdata}. Indices
#'   along the third dimension correspond to different resampling iterations.
#'   With multivariate outcomes, indices along the fourth dimension correspond
#'   to outcome-specific coefficients.} \item{method}{a list of
#'   \code{implementation}, \code{family}, \code{resampling} and
#'   \code{PFER_method} values used for the run.} \item{param}{a list of
#'   \code{K}, \code{pi_list}, \code{tau}, \code{n_cat}, \code{pk}, \code{n}
#'   (number of observations), \code{PFER_thr}, \code{FDP_thr} and \code{seed}
#'   values used for the run. The datasets \code{xdata} and \code{ydata} are
#'   also included if \code{output_data=TRUE}.} For all objects except those
#'   stored in \code{methods} or \code{params}, rows correspond to parameter
#'   values stored in the output \code{Lambda}.
#'
#' @keywords internal
SerialClustering <- function(xdata, Lambda, nc, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, tau = 0.5, seed = 1, n_cat = 3,
                             implementation = SparseHierarchicalClustering,
                             output_data = FALSE, verbose = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Defining resampling method (only subsampling is available as bootstrap would give distance of zero)
  resampling <- "subsampling"
  PFER_method <- "MB"
  PFER_thr <- Inf
  FDP_thr <- Inf

  # # Defining K if using complementary pairs (SS)
  # if (PFER_method == "SS") {
  #   K <- ceiling(K / 2) * 2
  #   tau <- 0.5
  # }

  # Initialising objects to be filled
  N <- N_block <- ncol(xdata)

  # Initialising array of co-membership proportions
  bigstab_obs <- array(0,
    dim = c(nrow(xdata), nrow(xdata), nrow(Lambda) * nrow(nc)),
    dimnames = list(rownames(xdata), rownames(xdata), NULL)
  )
  sampled_pairs <- matrix(0, nrow = nrow(xdata), ncol = nrow(xdata))
  rownames(sampled_pairs) <- colnames(sampled_pairs) <- rownames(xdata)

  # Initialising the array for contributing variables
  Beta <- array(0, dim = c(nrow(Lambda) * nrow(nc), ncol(xdata), K))
  rownames(Beta) <- paste0("s", seq(0, nrow(Beta) - 1))
  colnames(Beta) <- colnames(xdata)

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Computation of the selection proportions over Lambda
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  if (PFER_method == "MB") {
    for (k in 1:K) {
      # Subsampling observations
      s <- Resample(
        data = xdata, family = NULL, tau = tau, resampling = resampling, ...
      )
      Xsub <- xdata[s, ]

      # Applying clustering algorithm
      mybeta <- ClusteringAlgo(
        xdata = Xsub,
        Lambda = Lambda, nc = nc,
        implementation = implementation, ...
      )

      # Storing weights, used to define set of selected variables
      Beta[rownames(mybeta$weight), colnames(mybeta$weight), k] <- mybeta$weight

      # Storing co-membership status
      for (i in 1:dim(mybeta$comembership)[3]) {
        bigstab_obs[s, s, i] <- bigstab_obs[s, s, i] + mybeta$comembership[, , i]
      }

      # Storing sampled pairs
      sampled_pairs[s, s] <- sampled_pairs[s, s] + 1

      if (verbose) {
        utils::setTxtProgressBar(pb, k / K)
      }
    }

    # Computing the selection proportions
    bigstab_var <- matrix(NA, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab_var) <- colnames(Beta)
    rownames(bigstab_var) <- rownames(Beta)
    for (i in 1:nrow(Beta)) {
      for (j in 1:ncol(Beta)) {
        bigstab_var[i, j] <- sum(Beta[i, j, ] != 0) / K
      }
    }

    # Computing the co-membership proportions
    for (i in 1:dim(bigstab_obs)[3]) {
      bigstab_obs[, , i] <- bigstab_obs[, , i] / sampled_pairs
    }
    bigstab_obs[is.nan(bigstab_obs)] <- NA
  }

  # if (PFER_method == "SS") {
  #   for (k in 1:ceiling(K / 2)) {
  #     # Subsampling observations
  #     s <- Resample(
  #       data = xdata, family = NULL, tau = tau, resampling = resampling, ...
  #     )
  #     Xsub <- xdata[s, ]
  #
  #     # Applying clustering algorithm
  #     mybeta1 <- ClusteringAlgo(
  #       xdata = Xsub,
  #       Lambda = Lambda, nc = nc,
  #       implementation = implementation, ...
  #     )
  #
  #     # Complementary subset
  #     ns <- seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s]
  #     Xsub <- xdata[ns, ]
  #     mybeta2 <- ClusteringAlgo(
  #       xdata = Xsub,
  #       Lambda = Lambda, nc = nc,
  #       implementation = implementation, ...
  #     )
  #
  #     # Storing weights, used to define set of selected variables from first set
  #     Beta[rownames(mybeta1$weight), colnames(mybeta1$weight), k] <- mybeta1$weight
  #
  #     # Storing weights, used to define set of selected variables from complementary set
  #     Beta[rownames(mybeta2$weight), colnames(mybeta2$weight), ceiling(K / 2) + k] <- mybeta2$weight
  #
  #     # Storing co-membership status
  #     for (i in 1:dim(mybeta$comembership)[3]) {
  #       bigstab_obs[s, s, i] <- bigstab_obs[s, s, i] + mybeta1$comembership[, , i]
  #       bigstab_obs[ns, ns, i] <- bigstab_obs[ns, ns, i] + mybeta2$comembership[, , i]
  #     }
  #
  #     # Storing sampled pairs
  #     sampled_pairs[s, s] <- sampled_pairs[s, s] + 1
  #     sampled_pairs[ns, ns] <- sampled_pairs[ns, ns] + 1
  #
  #     if (verbose) {
  #       utils::setTxtProgressBar(pb, 2 * k / K)
  #     }
  #   }
  #
  #   # Computing the simultaneous selection proportions
  #   bigstab_var <- matrix(0, nrow = nrow(Beta), ncol = ncol(Beta))
  #   colnames(bigstab_var) <- colnames(Beta)
  #   rownames(bigstab_var) <- rownames(Beta)
  #   for (k in 1:ceiling(K / 2)) {
  #     A1 <- ifelse(Beta[, , k] != 0, yes = 1, no = 0)
  #     A2 <- ifelse(Beta[, , ceiling(K / 2) + k] != 0, yes = 1, no = 0)
  #     A <- A1 + A2
  #     A <- ifelse(A == 2, yes = 1, no = 0)
  #     bigstab_var <- bigstab_var + A
  #   }
  #   bigstab_var <- bigstab_var / ceiling(K / 2)
  #
  #   # Computing the co-membership proportions (cannot do simultaneous as different observations)
  #   for (i in 1:dim(bigstab_obs)[3]) {
  #     bigstab_obs[, , i] <- bigstab_obs[, , i] / sampled_pairs
  #   }
  #   bigstab_obs[is.nan(bigstab_obs)] <- NA
  # }

  if (verbose) {
    cat("\n")
  }

  # Computation of the stability score over Lambda and pi_list
  metrics <- StabilityMetrics(
    selprop = bigstab_obs, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
    Sequential_template = NULL, graph = TRUE,
    PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr
  )
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
    S = metrics$S,
    Lambda = cbind(rep(Lambda[, 1], each = nrow(nc))),
    nc = cbind(rep(nc, nrow(Lambda))),
    Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
    PFER = metrics$PFER, FDP = metrics$FDP,
    S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
    coprop = bigstab_obs,
    selprop = bigstab_var,
    Beta = Beta,
    methods = list(
      type = "clustering", implementation = myimplementation,
      resampling = myresampling, PFER_method = PFER_method
    ),
    params = list(
      K = K, pi_list = pi_list, tau = tau, n_cat = n_cat,
      pk = ncol(xdata), n = nrow(xdata),
      PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      seed = seed
    )
  )

  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata))
  }

  return(out)
}

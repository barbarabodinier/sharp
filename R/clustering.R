#' Consensus clustering
#'
#' Runs consensus clustering models with different combinations of parameters
#' controlling the number of clusters in the underlying algorithm and thresholds
#' in co-membership proportions. These two parameters are jointly calibrated by
#' maximising the stability score of the model (possibly under a constraint on
#' the expected number of falsely stably selected features). This function can
#' be used to identify stable clusters of observations sharing similar profiles.
#'
#' @inheritParams VariableSelection
#' @param xdata data matrix with observations as rows and variables as columns.
#' @param Lambda matrix of parameters controlling the number of clusters in the
#'   underlying algorithm specified in \code{implementation}. If \code{Lambda}
#'   is not provided, it is set to \code{seq(1, nrow(xdata))}.
#' @param implementation function to use for clustering. Possible functions
#'   include \code{\link{HierarchicalClustering}} (hierarchical clustering),
#'   \code{\link{KMeansClustering}} (k-means), \code{\link{GMMClustering}}
#'   (Gaussian Mixture Models), \code{\link{PAMClustering}} (Partioning Around
#'   Medoids), \code{\link{CLARAClustering}} (Clustering Large Applications) and
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
#'   \code{\link{HierarchicalClustering}}, \code{\link{KMeansClustering}},
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
#'   Lambda = 1:3, K = 5,
#'   verbose = FALSE
#' )
#' mymembership <- Clusters(stab)
#' }
#' \dontrun{
#' # Simulation of 15 observations belonging to 3 groups
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100,
#'   v_within = c(-1, -0.5), continuous = TRUE
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # Consensus clustering based on hierarchical clustering
#' stab <- Clustering(xdata = simul$data)
#' CalibrationPlot(stab, xlab = expression(italic(k)))
#' SelectionProportions(stab)
#' plot(Graph(Adjacency(stab), satellites = TRUE))
#' table(simul$theta, Clusters(stab))
#'
#' # Consensus clustering based on k-means clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = KMeansClustering
#' )
#' table(simul$theta, Clusters(stab))
#'
#' # Consensus clustering based on Gaussian Mixture Models
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = GMMClustering
#' )
#' table(simul$theta, Clusters(stab))
#'
#' # Consensus clustering based on PAM clustering
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = PAMClustering
#' )
#' table(simul$theta, Clusters(stab))
#'
#' # Consensus clustering based on CLARA algorithm
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = CLARAClustering,
#'   Lambda = 2:5,
#'   sampsize = 10
#' )
#' table(simul$theta, Clusters(stab))
#'
#' # Consensus clustering based on DBSCAN algorithm
#' stab <- Clustering(
#'   xdata = simul$data,
#'   implementation = DBSCANClustering,
#'   Lambda = seq(6, 10, by = 0.1)
#' )
#' CalibrationPlot(stab, xlab = "eps")
#' table(simul$theta, Clusters(stab))
#' }
#' @export
Clustering <- function(xdata, Lambda = NULL,
                       pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3,
                       implementation = HierarchicalClustering, scale = TRUE,
                       resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                       n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  # Visiting all possible numbers of clusters
  if (is.null(Lambda)) {
    Lambda <- cbind(seq(1, nrow(xdata)))
  }

  # Transposing of xdata for consistency with clustering literature (i.e. clustering of rows)
  xdata <- t(xdata)

  # Error and warning messages
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
    return(SerialGraphical(
      xdata = xdata, pk = pk, Lambda = Lambda, lambda_other_blocks = 0.1,
      pi_list = pi_list, K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
      implementation = implementation, scale = scale,
      resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
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

    # Removing graphical model specific outputs
    out$methods <- out$methods[-which(names(out$methods) %in% "start")]
    out$params <- out$params[-which(names(out$params) %in% c("lambda_other_blocks"))]
    out <- out[-which(names(out) %in% c("sign"))]

    # Updating n and pk (using transpose)
    n <- out$params$pk
    pk <- out$params$n
    out$params$pk <- pk
    out$params$n <- out$params$n
  }

  return(out)
}

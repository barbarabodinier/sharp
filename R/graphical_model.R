#' Stability selection graphical model
#'
#' Performs stability selection for graphical models. The underlying graphical
#' model (e.g. graphical LASSO) is run with different combinations of parameters
#' controlling the sparsity (e.g. penalty parameter) and thresholds in selection
#' proportions. These two hyper-parameters are jointly calibrated by
#' maximisation of the stability score.
#'
#' @inheritParams VariableSelection
#' @param xdata data matrix with observations as rows and variables as columns.
#'   For multi-block stability selection, the variables in data have to be
#'   ordered by group.
#' @param pk optional vector encoding the grouping structure. Only used for
#'   multi-block stability selection where \code{pk} indicates the number of
#'   variables in each group. If \code{pk=NULL}, single-block stability
#'   selection is performed.
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   If \code{Lambda=NULL} and \code{implementation=PenalisedGraphical},
#'   \code{\link{LambdaGridGraphical}} is used to define a relevant grid.
#'   \code{Lambda} can be provided as a vector or a matrix with
#'   \code{length(pk)} columns.
#' @param lambda_other_blocks optional vector of parameters controlling the
#'   level of sparsity in neighbour blocks for the multi-block procedure. To use
#'   jointly a specific set of parameters for each block,
#'   \code{lambda_other_blocks} must be set to \code{NULL} (not recommended).
#'   Only used for multi-block stability selection, i.e. if \code{length(pk)>1}.
#' @param implementation function to use for graphical modelling. If
#'   \code{implementation=PenalisedGraphical}, the algorithm implemented in
#'   \code{\link[glassoFast]{glassoFast}} is used for regularised estimation of
#'   a conditional independence graph. Alternatively, a user-defined function
#'   can be provided.
#' @param start character string indicating if the algorithm should be
#'   initialised at the estimated (inverse) covariance with previous penalty
#'   parameters (\code{start="warm"}) or not (\code{start="cold"}). Using
#'   \code{start="warm"} can speed-up the computations, but could lead to
#'   convergence issues (in particular with small \code{Lambda_cardinal}). Only
#'   used for \code{implementation=PenalisedGraphical} (see argument
#'   \code{"start"} in \code{\link[glassoFast]{glassoFast}}).
#' @param scale logical indicating if the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix should be used as input of
#'   \code{\link[glassoFast]{glassoFast}} if
#'   \code{implementation=PenalisedGraphical}. Otherwise, this argument must be
#'   used in the function provided in \code{implementation}.
#' @param lambda_max optional maximum value for the grid in penalty parameters.
#'   If \code{lambda_max=NULL}, the maximum value is set to the maximum
#'   covariance in absolute value. Only used if
#'   \code{implementation=PenalisedGraphical} and \code{Lambda=NULL}.
#' @param lambda_path_factor multiplicative factor used to define the minimum
#'   value in the grid.
#' @param max_density threshold on the density. The grid is defined such that
#'   the density of the estimated graph does not exceed max_density.
#'
#' @details In stability selection, a feature selection algorithm is fitted on
#'   \code{K} subsamples (or bootstrap samples) of the data with different
#'   parameters controlling the sparsity (\code{Lambda}). For a given (set of)
#'   sparsity parameter(s), the proportion out of the \code{K} models in which
#'   each feature is selected is calculated. Features with selection proportions
#'   above a threshold pi are considered stably selected. The stability
#'   selection model is controlled by the sparsity parameter(s) for the
#'   underlying algorithm, and the threshold in selection proportion:
#'
#'   \eqn{V_{\lambda, \pi} = \{ j: p_{\lambda}(j) \ge \pi \} }
#'
#'   These parameters can be calibrated by maximisation of a stability score
#'   (see \code{\link{ConsensusScore}} if \code{n_cat=NULL} or
#'   \code{\link{StabilityScore}} otherwise) calculated under the null
#'   hypothesis of equiprobability of selection.
#'
#'   It is strongly recommended to examine the calibration plot carefully to
#'   check that the grids of parameters \code{Lambda} and \code{pi_list} do not
#'   restrict the calibration to a region that would not include the global
#'   maximum (see \code{\link{CalibrationPlot}}). In particular, the grid
#'   \code{Lambda} may need to be extended when the maximum stability is
#'   observed on the left or right edges of the calibration heatmap. In some
#'   instances, multiple peaks of stability score can be observed. Simulation
#'   studies suggest that the peak corresponding to the largest number of
#'   selected features tend to give better selection performances. This is not
#'   necessarily the highest peak (which is automatically retained by the
#'   functions in this package). The user can decide to manually choose another
#'   peak.
#'
#'   To control the expected number of False Positives (Per Family Error Rate)
#'   in the results, a threshold \code{PFER_thr} can be specified. The
#'   optimisation problem is then constrained to sets of parameters that
#'   generate models with an upper-bound in PFER below \code{PFER_thr} (see
#'   Meinshausen and Bühlmann (2010) and Shah and Samworth (2013)).
#'
#'   Possible resampling procedures include defining (i) \code{K} subsamples of
#'   a proportion \code{tau} of the observations, (ii) \code{K} bootstrap samples
#'   with the full sample size (obtained with replacement), and (iii) \code{K/2}
#'   splits of the data in half for complementary pair stability selection (see
#'   arguments \code{resampling} and \code{cpss}). In complementary pair
#'   stability selection, a feature is considered selected at a given resampling
#'   iteration if it is selected in the two complementary subsamples.
#'
#'   To ensure reproducibility of the results, the starting number of the random
#'   number generator is set to \code{seed}.
#'
#'   For parallelisation, stability selection with different sets of parameters
#'   can be run on \code{n_cores} cores. Using \code{n_cores > 1} creates a
#'   \code{\link[future]{multisession}}. Alternatively,
#'   the function can be run manually with different \code{seed}s and all other
#'   parameters equal. The results can then be combined using
#'   \code{\link{Combine}}.
#'
#'   The generated network can be converted into
#'   \code{\link[igraph:igraph-package]{igraph}} object using
#'   \code{\link{Graph}}. The R package
#'   \code{\link[visNetwork:visDocumentation]{visNetwork}} can be used for
#'   interactive network visualisation (see examples in \code{\link{Graph}}).
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#'   \insertRef{stabilityselectionSS}{sharp}
#'
#'   \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{GraphicalLasso}{sharp}
#'
#' @return An object of class \code{graphical_model}. A list with: \item{S}{a
#'   matrix of the best stability scores for different (sets of) parameters
#'   controlling the level of sparsity in the underlying algorithm.}
#'   \item{Lambda}{a matrix of parameters controlling the level of sparsity in
#'   the underlying algorithm.} \item{Q}{a matrix of the average number of
#'   selected features by the underlying algorithm with different parameters
#'   controlling the level of sparsity.} \item{Q_s}{a matrix of the calibrated
#'   number of stably selected features with different parameters controlling
#'   the level of sparsity.} \item{P}{a matrix of calibrated thresholds in
#'   selection proportions for different parameters controlling the level of
#'   sparsity in the underlying algorithm.} \item{PFER}{a matrix of upper-bounds
#'   in PFER of calibrated stability selection models with different parameters
#'   controlling the level of sparsity.} \item{FDP}{a matrix of upper-bounds in
#'   FDP of calibrated stability selection models with different parameters
#'   controlling the level of sparsity.} \item{S_2d}{a matrix of stability
#'   scores obtained with different combinations of parameters. Columns
#'   correspond to different thresholds in selection proportions.}
#'   \item{PFER_2d}{a matrix of upper-bounds in FDP obtained with different
#'   combinations of parameters. Columns correspond to different thresholds in
#'   selection proportions. Only returned if \code{length(pk)=1}.}
#'   \item{FDP_2d}{a matrix of upper-bounds in PFER obtained with different
#'   combinations of parameters. Columns correspond to different thresholds in
#'   selection proportions. Only returned if \code{length(pk)=1}.}
#'   \item{selprop}{an array of selection proportions. Rows and columns
#'   correspond to nodes in the graph. Indices along the third dimension
#'   correspond to different parameters controlling the level of sparsity in the
#'   underlying algorithm.} \item{sign}{a matrix of signs of Pearson's
#'   correlations estimated from \code{xdata}.} \item{method}{a list with
#'   \code{type="graphical_model"} and values used for arguments
#'   \code{implementation}, \code{start}, \code{resampling}, \code{cpss} and
#'   \code{PFER_method}.} \item{params}{a list with values used for arguments
#'   \code{K}, \code{pi_list}, \code{tau}, \code{n_cat}, \code{pk}, \code{n}
#'   (number of observations in \code{xdata}), \code{PFER_thr}, \code{FDP_thr},
#'   \code{seed}, \code{lambda_other_blocks}, and \code{Sequential_template}.}
#'   The rows of \code{S}, \code{Lambda}, \code{Q}, \code{Q_s}, \code{P},
#'   \code{PFER}, \code{FDP}, \code{S_2d}, \code{PFER_2d} and \code{FDP_2d}, and
#'   indices along the third dimension of \code{selprop} are ordered in the same
#'   way and correspond to parameter values stored in \code{Lambda}. For
#'   multi-block inference, the columns of \code{S}, \code{Lambda}, \code{Q},
#'   \code{Q_s}, \code{P}, \code{PFER} and \code{FDP}, and indices along the
#'   third dimension of \code{S_2d} correspond to the different blocks.
#'
#' @family stability functions
#'
#' @seealso \code{\link{PenalisedGraphical}}, \code{\link{GraphicalAlgo}},
#'   \code{\link{LambdaGridGraphical}}, \code{\link{Resample}},
#'   \code{\link{StabilityScore}} \code{\link{Graph}}, \code{\link{Adjacency}},
#'
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = rep(7, 4))
#'
#' ## Single-block stability selection
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, nu_within = 0.1)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#' print(stab)
#'
#' # Calibration heatmap
#' CalibrationPlot(stab)
#'
#' # Visualisation of the results
#' summary(stab)
#' plot(stab)
#'
#' # Extraction of adjacency matrix or igraph object
#' Adjacency(stab)
#' Graph(stab)
#'
#'
#' ## Multi-block stability selection
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = c(10, 10))
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data, pk = c(10, 10), Lambda_cardinal = 10)
#' print(stab)
#'
#' # Calibration heatmap
#' # par(mfrow = c(1, 3))
#' CalibrationPlot(stab) # Producing three plots
#'
#' # Visualisation of the results
#' summary(stab)
#' plot(stab)
#'
#' # Multi-parameter stability selection (not recommended)
#' Lambda <- matrix(c(0.8, 0.6, 0.3, 0.5, 0.4, 0.3, 0.7, 0.5, 0.1), ncol = 3)
#' stab <- GraphicalModel(
#'   xdata = simul$data, pk = c(10, 10),
#'   Lambda = Lambda, lambda_other_blocks = NULL
#' )
#' stab$Lambda
#'
#'
#' ## Example with user-defined function: shrinkage estimation and selection
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, nu_within = 0.1)
#'
#' if (requireNamespace("corpcor", quietly = TRUE)) {
#'   # Writing user-defined algorithm in a portable function
#'   ShrinkageSelection <- function(xdata, Lambda, ...) {
#'     mypcor <- corpcor::pcor.shrink(xdata, verbose = FALSE)
#'     adjacency <- array(NA, dim = c(nrow(mypcor), ncol(mypcor), nrow(Lambda)))
#'     for (k in 1:nrow(Lambda)) {
#'       A <- ifelse(abs(mypcor) >= Lambda[k, 1], yes = 1, no = 0)
#'       diag(A) <- 0
#'       adjacency[, , k] <- A
#'     }
#'     return(list(adjacency = adjacency))
#'   }
#'
#'   # Running the algorithm without stability
#'   myglasso <- GraphicalAlgo(
#'     xdata = simul$data,
#'     Lambda = matrix(c(0.05, 0.1), ncol = 1), implementation = ShrinkageSelection
#'   )
#'
#'   # Stability selection using shrinkage estimation and selection
#'   stab <- GraphicalModel(
#'     xdata = simul$data, Lambda = matrix(c(0.01, 0.05, 0.1), ncol = 1),
#'     implementation = ShrinkageSelection
#'   )
#'   CalibrationPlot(stab)
#'   stable_adjacency <- Adjacency(stab)
#' }
#'
#' par(oldpar)
#' }
#' @export
GraphicalModel <- function(xdata, pk = NULL, Lambda = NULL, lambda_other_blocks = 0.1,
                           pi_list = seq(0.01, 0.99, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = NULL,
                           implementation = PenalisedGraphical, start = "warm", scale = TRUE,
                           resampling = "subsampling", cpss = FALSE,
                           PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                           Lambda_cardinal = 50, lambda_max = NULL, lambda_path_factor = 0.001, max_density = 0.5,
                           n_cores = 1, output_data = FALSE, verbose = TRUE, beep = NULL, ...) {
  # Definition of the type of approach (single or multi-block)
  if (is.null(pk)) {
    pk <- ncol(xdata)
  }
  if (length(pk) > 1) {
    calibration <- "multi-block"
  } else {
    calibration <- "single-block"
  }

  # Error and warning messages
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

  # Launching stability selection and calibration
  if (is.null(Lambda)) {
    # Defining a broad grid of lambda values
    Lambda <- LambdaGridGraphical(
      xdata = xdata, pk = pk, lambda_other_blocks = lambda_other_blocks, tau = tau,
      implementation = implementation, start = "cold", scale = scale,
      resampling = resampling, cpss = cpss,
      PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      lambda_max = lambda_max, lambda_path_factor = lambda_path_factor, max_density = max_density,
      Lambda_cardinal = Lambda_cardinal, ...
    )
  }

  # Stability selection and score
  if (n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
    mypar <- future.apply::future_lapply(X = 1:n_cores, future.seed = TRUE, FUN = function(k) {
      return(SerialGraphical(
        xdata = xdata, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
        pi_list = pi_list, K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
        implementation = implementation, start = start, scale = scale,
        resampling = resampling, cpss = cpss, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
        output_data = output_data, verbose = FALSE, ...
      ))
    })
    future::plan(future::sequential)

    # Combining the outputs from parallel iterations
    out <- mypar[[1]]
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[i]]))
    }
  } else {
    out <- SerialGraphical(
      xdata = xdata, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
      pi_list = pi_list, K = K, tau = tau, seed = seed, n_cat = n_cat,
      implementation = implementation, start = start, scale = scale,
      resampling = resampling, cpss = cpss, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
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
  }

  # Defining the class
  class(out) <- "graphical_model"

  # Making beep
  if (!is.null(beep)) {
    beepr::beep(sound = beep)
  }

  return(out)
}


#' Stability selection graphical model (internal)
#'
#' Runs stability selection graphical models with different combinations of
#' parameters controlling the sparsity of the underlying selection algorithm
#' (e.g. penalty parameter for regularised models) and thresholds in selection
#' proportions. These two parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features). This function uses a serial
#' implementation and requires the grid of parameters controlling the underlying
#' algorithm as input (for internal use only).
#'
#' @inheritParams GraphicalModel
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   If \code{implementation="glassoFast"}, \code{Lambda} contains penalty
#'   parameters.
#'
#' @return A list with: \item{S}{a matrix of the best stability scores for
#'   different (sets of) parameters controlling the level of sparsity in the
#'   underlying algorithm.} \item{Lambda}{a matrix of parameters controlling the
#'   level of sparsity in the underlying algorithm.} \item{Q}{a matrix of the
#'   average number of selected features by the underlying algorithm with
#'   different parameters controlling the level of sparsity.} \item{Q_s}{a
#'   matrix of the calibrated number of stably selected features with different
#'   parameters controlling the level of sparsity.} \item{P}{a matrix of
#'   calibrated thresholds in selection proportions for different parameters
#'   controlling the level of sparsity in the underlying algorithm.}
#'   \item{PFER}{a matrix of upper-bounds in PFER of calibrated stability
#'   selection models with different parameters controlling the level of
#'   sparsity.} \item{FDP}{a matrix of upper-bounds in FDP of calibrated
#'   stability selection models with different parameters controlling the level
#'   of sparsity.} \item{S_2d}{a matrix of stability scores obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in selection proportions.} \item{PFER_2d}{a matrix of
#'   upper-bounds in FDP obtained with different combinations of parameters.
#'   Columns correspond to different thresholds in selection proportions. Only
#'   returned if \code{length(pk)=1}.} \item{FDP_2d}{a matrix of upper-bounds in
#'   PFER obtained with different combinations of parameters. Columns correspond
#'   to different thresholds in selection proportions. Only returned if
#'   \code{length(pk)=1}.} \item{selprop}{an array of selection proportions.
#'   Rows and columns correspond to nodes in the graph. Indices along the third
#'   dimension correspond to different parameters controlling the level of
#'   sparsity in the underlying algorithm.} \item{sign}{a matrix of signs of
#'   Pearson's correlations estimated from \code{xdata}.} \item{method}{a list
#'   with \code{type="graphical_model"} and values used for arguments
#'   \code{implementation}, \code{start}, \code{resampling}, \code{cpss} and
#'   \code{PFER_method}.} \item{params}{a list with values used for arguments
#'   \code{K}, \code{pi_list}, \code{tau}, \code{n_cat}, \code{pk}, \code{n}
#'   (number of observations in \code{xdata}), \code{PFER_thr}, \code{FDP_thr},
#'   \code{seed}, \code{lambda_other_blocks}, and \code{Sequential_template}.}
#'   The rows of \code{S}, \code{Lambda}, \code{Q}, \code{Q_s}, \code{P},
#'   \code{PFER}, \code{FDP}, \code{S_2d}, \code{PFER_2d} and \code{FDP_2d}, and
#'   indices along the third dimension of \code{selprop} are ordered in the same
#'   way and correspond to parameter values stored in \code{Lambda}. For
#'   multi-block inference, the columns of \code{S}, \code{Lambda}, \code{Q},
#'   \code{Q_s}, \code{P}, \code{PFER} and \code{FDP}, and indices along the
#'   third dimension of \code{S_2d} correspond to the different blocks.
#'
#' @keywords internal
SerialGraphical <- function(xdata, pk = NULL, Lambda, lambda_other_blocks = 0.1,
                            pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = n_cat,
                            implementation = PenalisedGraphical, start = "cold", scale = TRUE,
                            resampling = "subsampling", cpss = FALSE,
                            PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                            output_data = FALSE, verbose = TRUE, ...) {
  # Marginal correlation to get sign of the relationship
  mycor_for_sign <- stats::cor(xdata)

  # Using complementary pairs for SS
  if (PFER_method == "SS") {
    cpss <- TRUE
  }

  # Defining K if using complementary pairs
  if (cpss) {
    K <- ceiling(K / 2) * 2
    tau <- 0.5
  }

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  nblocks <- length(pk) * (length(pk) + 1) / 2
  bigblocks_vect <- factor(bigblocks[upper.tri(bigblocks)], levels = 1:nblocks)
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- levels(bigblocks_vect)
  names(N_blocks) <- blocks

  # Preparing the PFER and FDP thresholds
  if (length(PFER_thr) == 1) {
    PFER_thr_blocks <- ceiling(prop.table(N_blocks) * PFER_thr)
    PFER_thr_blocks[which(is.na(PFER_thr_blocks))] <- Inf
  } else {
    if (length(PFER_thr) == nblocks) {
      PFER_thr_blocks <- PFER_thr
    }
  }
  if (length(FDP_thr) == 1) {
    FDP_thr_blocks <- rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr) == nblocks) {
      FDP_thr_blocks <- FDP_thr
    }
  }

  # Re-formatting Lambda
  if (is.vector(Lambda)) {
    grid <- BlockLambdaGrid(Lambda = Lambda, lambda_other_blocks = lambda_other_blocks)
    Lambda <- grid$Lambda
    Sequential_template <- grid$Sequential_template
  } else {
    grid <- BlockLambdaGrid(Lambda = Lambda, lambda_other_blocks = lambda_other_blocks)
    Lambda <- grid$Lambda
    Sequential_template <- grid$Sequential_template
  }

  # Initialising array of selection proportions
  bigstab <- array(0,
    dim = c(ncol(xdata), ncol(xdata), nrow(Lambda)),
    dimnames = list(colnames(xdata), colnames(xdata), NULL)
  )

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Initialisation of the run
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }

  # Using MB formula of the PFER
  if (!cpss) {
    for (i in 1:K) {
      # Resampling of the data
      s <- Resample(data = xdata, family = NULL, tau = tau, resampling = resampling, ...)
      xdata_sub <- xdata[s, , drop = FALSE]

      # Estimation of the networks for different penalties
      A <- GraphicalAlgo(
        xdata = xdata_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
        scale = scale, implementation = implementation, start = start, ...
      )

      # Computing the selection counts
      for (k in 1:dim(A)[3]) {
        bigstab[, , k] <- bigstab[, , k] + A[, , k]
      }

      if (verbose) {
        utils::setTxtProgressBar(pb, i / K)
      }
    }
    # Getting selection proportions from selection counts
    for (k in 1:dim(bigstab)[3]) {
      bigstab[, , k] <- bigstab[, , k] / K
      diag(bigstab[, , k]) <- 0
    }
  } else {
    for (i in 1:ceiling(K / 2)) {
      # Sample 1
      s <- Resample(data = xdata, family = NULL, tau = tau, resampling = resampling, ...)
      xdata_sub <- xdata[s, , drop = FALSE]

      # Estimation of the networks for different penalties
      A1 <- GraphicalAlgo(
        xdata = xdata_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
        scale = scale, implementation = implementation, start = start, ...
      )

      # Sample 2: everything not in sample 1
      xdata_sub <- xdata[-s, , drop = FALSE]

      # Estimation of the networks for different penalties
      A2 <- GraphicalAlgo(
        xdata = xdata_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
        scale = scale, implementation = implementation, start = start, ...
      )

      # Computing the simultaneous selection counts
      for (k in 1:dim(A1)[3]) {
        A <- ifelse((A1[, , k] + A2[, , k]) == 2, yes = 1, no = 0)
        bigstab[, , k] <- bigstab[, , k] + A
      }

      if (verbose) {
        utils::setTxtProgressBar(pb, i / ceiling(K / 2))
      }
    }
    # Getting selection proportions from selection counts
    for (k in 1:dim(bigstab)[3]) {
      bigstab[, , k] <- bigstab[, , k] / ceiling(K / 2)
      diag(bigstab[, , k]) <- 0
    }
  }

  # Computation of the stability score
  metrics <- StabilityMetrics(
    selprop = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
    Sequential_template = Sequential_template, graph = TRUE,
    PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
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
  if (nblocks == 1) {
    out <- list(
      S = metrics$S, Lambda = Lambda,
      Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
      PFER = metrics$PFER, FDP = metrics$FDP,
      S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
      selprop = bigstab, sign = sign(mycor_for_sign),
      methods = list(
        type = "graphical_model", implementation = myimplementation, start = start,
        resampling = myresampling, cpss = cpss, PFER_method = PFER_method
      ),
      params = list(
        K = K, pi_list = pi_list, tau = tau, n_cat = n_cat,
        pk = pk, n = nrow(xdata),
        PFER_thr = PFER_thr, FDP_thr = FDP_thr, seed = seed,
        lambda_other_blocks = lambda_other_blocks, Sequential_template = Sequential_template
      )
    )
    if (output_data) {
      out$params <- c(out$params, list(xdata = xdata))
    }
  } else {
    out <- list(
      S = metrics$S, Lambda = Lambda,
      Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
      PFER = metrics$PFER, FDP = metrics$FDP,
      S_2d = metrics$S_2d,
      selprop = bigstab, sign = sign(mycor_for_sign),
      methods = list(
        type = "graphical_model", implementation = myimplementation, start = start,
        resampling = myresampling, cpss = cpss, PFER_method = PFER_method
      ),
      params = list(
        K = K, pi_list = pi_list, tau = tau, n_cat = n_cat,
        pk = pk, n = nrow(xdata),
        PFER_thr = PFER_thr, FDP_thr = FDP_thr, seed = seed,
        lambda_other_blocks = lambda_other_blocks, Sequential_template = Sequential_template
      )
    )
    if (output_data) {
      out$params <- c(out$params, list(xdata = xdata))
    }
  }

  # Defining the class
  class(out) <- "graphical_model"

  return(out)
}

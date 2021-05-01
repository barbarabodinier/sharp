#' Stability selection graphical model
#'
#' Runs stability selection graphical models with different combinations of
#' parameters controlling the sparsity of the underlying selection algorithm
#' (e.g. penalty parameter for regularised models) and thresholds in selection
#' proportions. These two parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features).
#'
#' @inheritParams VariableSelection
#' @param data matrix with observations as rows and variables as columns. For
#'   multi-block stability selection, the variables in data have to be ordered
#'   by group.
#' @param pk optional vector encoding the grouping structure. Only used for
#'   multi-block stability selection where \code{pk} indicates the number of
#'   variables in each group. If \code{pk=NULL}, single-block stability
#'   selection is performed.
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   If \code{implementation="glassoFast"}, \code{Lambda} contains penalty
#'   parameters. If \code{Lambda=NULL}, \code{\link{LambdaGridGraphical}} is
#'   used to define a relevant grid. \code{Lambda} can be provided as a vector
#'   or a matrix with \code{length(pk)} columns. If \code{implementation} is not
#'   set to \code{"glassoFast"}, \code{Lambda} must be provided.
#' @param lambda_other_blocks optional vector of parameters controlling the
#'   level of sparsity in neighbour blocks for the multi-block procedure. To use
#'   jointly a specific set of parameters for each block,
#'   \code{lambda_other_blocks} must be set to \code{NULL} (not recommended).
#'   Only used for multi-block stability selection, i.e. if \code{length(pk)>1}.
#' @param implementation character string indicating the name of the function to
#'   use for graphical modelling. With \code{implementation="glassoFast"},
#'   \code{\link[glassoFast]{glassoFast}} is used for regularised estimation of
#'   a conditional independence graph. Alternatively, a function with arguments
#'   \code{x}, \code{lambda}, \code{scale} and \code{...}, and returning a
#'   binary and symmetric matrix for which diagonal elements are equal to zero
#'   can be used (more details in \code{\link{GraphicalAlgo}}).
#' @param start character string indicating if the algorithm should be
#'   initialised at the estimated (inverse) covariance with previous penalty
#'   parameters (\code{start="warm"}) or not (\code{start="cold"}). Using
#'   \code{start="warm"} can speed-up the computations. Only used for
#'   \code{implementation="glassoFast"} (see argument \code{"start"} in
#'   \code{\link[glassoFast]{glassoFast}}).
#' @param scale logical indicating if the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix should be used as input for the
#'   graphical LASSO if \code{implementation="glassoFast"}. Otherwise, this
#'   argument must be used in the function provided in \code{implementation}.
#' @param lambda_max optional maximum value for the grid in penalty parameters.
#'   If \code{lambda_max=NULL}, the maximum value is set to the maximum
#'   covariance in absolute value. Only used if
#'   \code{implementation="glassoFast"} and \code{Lambda=NULL}.
#' @param lambda_path_factor multiplicative factor used to define the minimum
#'   value in the grid.
#' @param max_density threshold on the density. The grid is defined such that
#'   the density of the estimated graph does not exceed max_density.
#'
#' @details To ensure reproducibility of the results, the state of the random
#'   number generator is fixed to \code{seed}. For parallelisation of the code,
#'   stability selection results produced with different \code{seed}s and all
#'   other parameters equal can be combined (more details in
#'   \code{\link{Combine}}).
#'
#' @return A list with: \item{S}{a matrix of the best (block-specific) stability
#'   scores for different (sets of) parameters controlling the level of sparsity
#'   in the underlying algorithm. } \item{Lambda}{a matrix of (block-specific)
#'   parameters controlling the level of sparsity. } \item{Q}{a matrix of
#'   average numbers of (block-specific) edges selected by the underlying
#'   algorihm for different (sets of) parameters controlling the level of
#'   sparsity.} \item{Q_s}{a matrix of calibrated numbers of (block-specific)
#'   stable edges for different (sets of) parameters controlling the level of
#'   sparsity in the underlying algorithm. } \item{P}{a matrix of calibrated
#'   (block-specific) thresholds in selection proportions for different (sets
#'   of) parameters controlling the level of sparsity in the underlying
#'   algorithm. } \item{PFER}{a matrix of the (block-specific) upper-bounds in
#'   PFER of calibrated stability selection models with different (sets of)
#'   parameters controlling the level of sparsity in the underlying algorithm.}
#'   \item{FDP}{a matrix of the (block-specific) upper-bounds in FDP of
#'   calibrated stability selection models with different (sets of) parameters
#'   controlling the level of sparsity in the underlying algorithm.}
#'   \item{S_2d}{an array of (block-specific) stability scores obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in selection proportions. In multi-block stability selection,
#'   indices along the third dimension correspond to different blocks.}
#'   \item{PFER_2d}{an array of computed upper-bounds of PFER obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in selection proportions. Only available for single-block
#'   stability selection.} \item{FDP_2d}{an array of computed upper-bounds of
#'   FDP obtained with different combinations of parameters. Columns correspond
#'   to different thresholds in selection proportions. Only available for
#'   single-block stability selection.} \item{selprop}{an array of selection
#'   proportions. Rows and columns correspond to nodes in the graph. Indices
#'   along the third dimension correspond to different (sets of) parameters
#'   controlling the level of sparsity in the underlying algorithm.}
#'   \item{sign}{a matrix of signs of Pearson's correlations estimated from
#'   \code{data}.} \item{method}{a list with \code{implementation},
#'   \code{start}, \code{resampling} and \code{PFER_method} values used for the
#'   run.} \item{param}{a list with values of other objects used for the run.}
#'   For all objects except \code{selprop},
#'   \code{sign} and those stored in \code{methods} or \code{params}, rows
#'   correspond to parameter values stored in the output \code{Lambda}. In
#'   multi-block stability selection, columns of these same objects except
#'   correspond to different blocks.
#'
#' @family stability selection functions
#' @seealso \code{\link{LambdaGridGraphical}}, \code{\link{Resample}},
#'   \code{\link{GraphicalAlgo}}
#'
#' @examples
#' \dontshow{
#' # Single-block stability selection
#' set.seed(1)
#' simul <- SimulateGraphical(n = 50, pk = 10, nu = 0.1)
#' stab <- GraphicalModel(data = simul$data, K = 5, verbose = FALSE)
#' A <- Adjacency(stab)
#' mygraph <- Graph(A)
#' perf <- SelectionPerformance(theta = A, theta_star = simul$theta)
#' perfgraph <- SelectionPerformanceGraph(
#'   theta = Adjacency(stab),
#'   theta_star = simul$theta, plot = TRUE
#' )
#'
#' # Multi-block stability selection
#' set.seed(1)
#' pk <- c(10, 10)
#' simul <- SimulateGraphical(n = 50, pk = pk)
#' stab <- GraphicalModel(data = simul$data, pk = pk, Lambda_cardinal = 10, K = 5, verbose = FALSE)
#' A <- Adjacency(stab)
#' mygraph <- Graph(A)
#' perf <- SelectionPerformance(theta = A, theta_star = simul$theta, pk = pk)
#' }
#' \dontrun{
#'
#' # Single-block stability selection
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, nu = 0.1)
#' stab <- GraphicalModel(data = simul$data)
#' plot(Graph(stab))
#'
#' # Multi-block stability selection
#' set.seed(1)
#' simul <- SimulateGraphical(pk = c(10, 10))
#' stab <- GraphicalModel(data = simul$data, pk = c(10, 10), Lambda_cardinal = 10)
#' stab$Lambda # sets of penalty parameters used jointly
#'
#' # Multi-parameter stability selection (not recommended)
#' Lambda <- matrix(c(0.8, 0.6, 0.3, 0.5, 0.4, 0.3, 0.7, 0.5, 0.1), ncol = 3)
#' stab <- GraphicalModel(
#'   data = simul$data, pk = c(10, 10),
#'   Lambda = Lambda, lambda_other_blocks = NULL
#' )
#' stab$Lambda
#' }
#' @export
GraphicalModel <- function(data, pk = NULL, Lambda = NULL, lambda_other_blocks = 0.1,
                           pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3,
                           implementation = "glassoFast", start = "warm", scale = TRUE,
                           resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                           Lambda_cardinal = 50, lambda_max = NULL, lambda_path_factor = 0.001, max_density = 0.5,
                           n_cores = 1, verbose = TRUE, ...) {
  # Definition of the type of approach (single or multi-block)
  if (is.null(pk)) {
    pk <- ncol(data)
  }
  if (length(pk) > 1) {
    calibration <- "multi-block"
  } else {
    calibration <- "single-block"
  }
  if (verbose) {
    print(paste("Starting", calibration, "calibration..."))
  }

  # Error and warning messages
  bigblocks <- bigblocks_vect <- blocks <- N_blocks <- nblocks <- PFER_thr_blocks <- FDP_thr_blocks <- NULL
  CheckInputGraphical(
    data = data, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
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
    if (verbose) {
      print("Defining the grid of lambda values...")
    }
    Lambda <- LambdaGridGraphical(
      data = data, pk = pk, lambda_other_blocks = lambda_other_blocks, tau = tau,
      implementation = implementation, start = "cold", scale = scale,
      resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      lambda_max = lambda_max, lambda_path_factor = lambda_path_factor, max_density = max_density,
      Lambda_cardinal = Lambda_cardinal, ...
    )
  }

  # Check if parallelisation is possible (forking)
  if (.Platform$OS.type != "unix") {
    if (n_cores > 1) {
      warning("Invalid input for argument 'n_cores'. Parallelisation relies on forking, it is not available on Windows.")
    }
    n_cores <- 1
  }

  # Stability selection and score
  mypar <- parallel::mclapply(X = 1:n_cores, FUN = function(k) {
    return(SerialGraphical(
      data = data, pk = pk, Lambda = Lambda, lambda_other_blocks = lambda_other_blocks,
      pi_list = pi_list, K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
      implementation = implementation, start = start, scale = scale,
      resampling = resampling, PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      verbose = verbose, ...
    ))
  })

  # Combining the outputs from parallel iterations
  out <- mypar[[1]]
  if (n_cores > 1) {
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[2]], graph = TRUE))
    }
  }

  if (verbose) {
    cat("\n")
    print("Visited Q:")
    if (nrow(out$Q) > 15) {
      print(utils::head(out$Q))
      print("[...]")
      print(utils::tail(out$Q))
    } else {
      print(out$Q)
    }
  }

  return(out)
}


#' Stability selection graphical model (internal)
#'
#' Runs stability selection graphical models with
#' different combinations of parameters controlling the sparsity
#' of the underlying selection algorithm (e.g. penalty parameter for
#' regularised models) and thresholds in selection proportions.
#' These two parameters are jointly calibrated by maximising
#' the stability score of the model (possibly under a constraint
#' on the expected number of falsely stably selected features).
#' This function uses a serial implementation and requires
#' the grid of parameters controlling the underlying algorithm
#' as input (for internal use only).
#'
#' @param data matrix with observations as rows and variables as columns.
#' @param pk vector encoding the grouping structure.
#' Only used for multi-block stability selection.
#' For this, the variables in data have to be ordered
#' by group and argument "pk" has to be a vector
#' indicating the number of variables
#' in each of the groups (see example below).
#' If pk=NULL, single-block stability selection is performed.
#' @param Lambda matrix of parameters controlling the underlying
#' feature selection algorithm specified in "implementation".
#' With implementation="glassoFast", these are penalty parameters
#' controlling the regularised model.
#' If Lambda=NULL, \code{\link{LambdaGridGraphical}} is used to define
#' a relevant grid.
#' For multi-block calibration (i.e. when argument "pk" is a vector),
#' Lambda can be a vector, to use the procedure from Equation (5) (recommended),
#' or a matrix with as many columns as there are entries in "pk",
#' to use the procedure from Equation (4) (see details and examples below).
#' @param lambda_other_blocks optional vector of (penalty) parameters
#' to use for other blocks
#' in the iterative multi-block procedure (see example below).
#' Only used for multi-block graphical models, i.e. when pk is a vector.
#' @param pi_list grid of values for the threshold in selection proportion.
#' With n_cat=3, these values must be between 0.5 and 1.
#' With n_cat=2, these values must be between 0 and 1.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param seed value of the seed to use to ensure reproducibility.
#' @param n_cat number of categories used to compute the stability score.
#' Possible values are 2 or 3.
#' @param implementation name of the function to use for definition of the grid
#' of lambda values. With implementation="glassoFast", the function \code{\link[glassoFast]{glassoFast}}
#' is called and iteratively applied on possible penalty values until the constraint are verified,
#' i.e. that the expected density is below the value given in "max_density",
#' that the expected PFER is below the value given in "PFER_thr"
#' or that the expected PFER is below the number of selected edges if "FDP_thr" is not set to Inf.
#' Alternatively, this argument can be a character string indicating the name of a function.
#' The function provided must use arguments called "x", "lambda" and "scale"
#' and return a binary and symmetric adjacency matrix (see example below).
#' @param start character string indicating if the algorithm should be
#' initialised at the estimated (inverse) covariance with previous
#' penalty parameters (start="warm") or not (start="cold").
#' Using start="warm" can speed-up the computations.
#' Only used for implementation="glassoFast" (see argument "start"
#' in \code{\link[glassoFast]{glassoFast}}).
#' @param scale logical indicating if the correlation (if scale=TRUE)
#' or covariance (if scale=FALSE) matrix should be used as input
#' for the graphical LASSO. If implementation is not set to "glassoFast",
#' this argument must be used as input of the function provided instead.
#' @param resampling resampling approach. Possible values are: "subsampling"
#' for sampling without replacement of a proportion tau of the observations, or
#' "bootstrap" for sampling with replacement generating a resampled dataset with
#' as many observations as in the full sample. Alternatively, this argument can be
#' a character string indicating the name of a function to use for resampling.
#' This function must use arguments called "data" and "tau" and return
#' IDs of observations to be included in the resampled dataset
#' (see example in \code{\link{Resample}}).
#' @param PFER_method method used to compute the expected number of False Positives,
#' (or Per Family Error Rate, PFER). With PFER_method="MB", the method
#' proposed by Meinshausen and Buhlmann (2010) is used. With PFER_method="SS",
#' the method proposed by Shah and Samworth (2013) under the assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error control.
#' With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is used.
#' The grid is defined such that the estimated graph does not generate an upper-bound in
#' PFER above PFER_thr.
#' @param FDP_thr threshold in the expected proportion of falsely selected edges
#' (or False Discovery Proportion, FDP)
#' for constrained calibration by error control.
#' With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is used.
#' If FDP_thr is not infinite, the grid is defined such that the estimated graph
#' does not generate an upper-bound in PFER above the number of node pairs.
#' @param verbose logical indicating if a message with minimum and maximum
#' numbers of selected variables on one instance of resampled data should be printed.
#' @param ... additional parameters passed to the functions provided in
#' "implementation" or "resampling".
#'
#' @return A list with:
#' \item{S}{a matrix of
#' the best (block-specific) stability scores
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{Lambda}{a matrix of
#' (block-specific) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to sets of penalty parameters and
#' columns correspond to different blocks.}
#' \item{Q}{a matrix of
#' average numbers of (block-specific) edges
#' selected by the underlying algorihm
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{Q_s}{a matrix of
#' calibrated numbers of (block-specific) stable edges
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{P}{a matrix of
#' calibrated (block-specific) thresholds in selection proportions
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{PFER}{a matrix of
#' computed (block-specific) upper-bounds in PFER of
#' calibrated graphs
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{FDP}{a matrix of
#' computed (block-specific) upper-bounds in FDP of
#' calibrated stability selection models
#' for different (sets of) penalty parameters.
#' In multi-block stability selection,
#' rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda") and
#' columns correspond to different blocks.}
#' \item{S_2d}{an array of
#' (block-specific) stability scores obtained
#' with different combinations of parameters.
#' Rows correspond to different (sets of) penalty parameters and
#' columns correspond to different thresholds in selection proportions.
#' In multi-block stability selection,
#' indices along the third dimension
#' correspond to different blocks.}
#' \item{PFER_2d}{an array of
#' computed upper-bounds of PFER obtained
#' with different combinations of parameters.
#' Rows correspond to different penalty parameters and
#' columns correspond to different thresholds in selection proportions.
#' Only available for single-block stability selection.}
#' \item{FDP_2d}{an array of
#' computed upper-bounds of FDP obtained
#' with different combinations of parameters.
#' Rows correspond to different penalty parameters and
#' columns correspond to different thresholds in selection proportions.
#' Only available for single-block stability selection.}
#' \item{selprop}{an array of selection proportions.
#' Rows and columns correspond to nodes in the graph.
#' Indices along the third dimension correspond to
#' different (sets of) penalty parameters.}
#' \item{sign}{a matrix of signs of Pearson's correlations.}
#' \item{method}{a list with input values for the arguments
#' "implementation", "start", "resampling" and "PFER_method".}
#' \item{param}{a list with input values for the arguments
#' "K", "pi_list", "tau", "n_cat", "pk", "PFER_thr", "FDP_thr",
#' "seed", "lambda_other_blocks" and "data".
#' The object "Sequential_template" is also returned
#' (for internal use).}
#'
#' @keywords internal
SerialGraphical <- function(data, pk = NULL, Lambda, lambda_other_blocks = 0.1,
                            pi_list = seq(0.6, 0.9, by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = n_cat,
                            implementation = "glassoFast", start = "cold", scale = TRUE,
                            resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                            verbose = TRUE, ...) {
  # Marginal correlation to get sign of the relationship
  mycor_for_sign <- stats::cor(data)

  # Defining K if using complementary pairs (SS)
  if (PFER_method == "SS") {
    K <- ceiling(K / 2) * 2
    tau <- 0.5
  }

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

  # Preparing the PFER and FDP thresholds
  if (length(PFER_thr) == 1) {
    PFER_thr_blocks <- ceiling(prop.table(N_blocks) * PFER_thr)
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

  # Showing the grid of (block-specific) lambda values
  if (verbose) {
    print("Grid of lambda values:")
    if (ncol(Lambda) == 1) {
      print(as.vector(Lambda))
    } else {
      print(Lambda)
    }
  }

  # Initialising array of selection proportions
  bigstab <- array(0,
    dim = c(ncol(data), ncol(data), nrow(Lambda)),
    dimnames = list(colnames(data), colnames(data), NULL)
  )

  # Printing message
  if (verbose) {
    if (all(!is.infinite(PFER_thr_blocks))) {
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))) {
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Initialisation of the run
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }

  # Using MB formula of the PFER
  if (PFER_method == "MB") {
    for (i in 1:K) {
      # Subsampling of the data
      s <- Resample(data = data, family = NULL, tau = tau, resampling = resampling, ...)
      data_sub <- data[s, ]

      # Estimation of the networks for different penalties
      A <- GraphicalAlgo(
        x = data_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
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
  }

  # Using complementary pairs and SS formula of the PFER
  if (PFER_method == "SS") {
    for (i in 1:ceiling(K / 2)) {
      # Sample 1
      s <- Resample(data = data, family = NULL, tau = tau, resampling = resampling, ...)
      data_sub <- data[s, ]

      # Estimation of the networks for different penalties
      A1 <- GraphicalAlgo(
        x = data_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
        scale = scale, implementation = implementation, start = start, ...
      )

      # Sample 2: everything not in sample 1
      data_sub <- data[-s, ]

      # Estimation of the networks for different penalties
      A2 <- GraphicalAlgo(
        x = data_sub, pk = pk, Lambda = Lambda, Sequential_template = Sequential_template,
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
  if (K > 2) {
    metrics <- StabilityMetrics(
      bigstab = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = TRUE,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
    if (verbose) {
      utils::setTxtProgressBar(pb, 1)
      cat("\n")
    }
  } else {
    # Initialising objects to be filled
    Q <- matrix(NA, nrow = nrow(Lambda), ncol = nblocks)
    for (k in 1:nrow(Lambda)) {
      # Extracting corresponding selection proportions
      stab_iter <- bigstab[, , k]

      # Getting number of selected variables per block
      for (block_id in 1:nblocks) {
        stab_iter_block <- stab_iter[(bigblocks == block_id) & (upper.tri(bigblocks))] # selection proportions in the block
        q_block <- round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
        Q[k, block_id] <- q_block
      }
    }
  }

  # Preparing outputs
  if (K > 2) {
    if (nblocks == 1) {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab, sign = sign(mycor_for_sign),
        methods = list(implementation = implementation, start = start, resampling = resampling, PFER_method = PFER_method),
        params = list(
          K = K, pi_list = pi_list, tau = tau, n_cat = n_cat, pk = pk, PFER_thr = PFER_thr, FDP_thr = FDP_thr, seed = seed,
          lambda_other_blocks = lambda_other_blocks, Sequential_template = Sequential_template, data = data
        )
      ))
    } else {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d,
        selprop = bigstab, sign = sign(mycor_for_sign),
        methods = list(implementation = implementation, start = start, resampling = resampling, PFER_method = PFER_method),
        params = list(
          K = K, pi_list = pi_list, tau = tau, n_cat = n_cat, pk = pk, PFER_thr = PFER_thr, FDP_thr = FDP_thr, seed = seed,
          lambda_other_blocks = lambda_other_blocks, Sequential_template = Sequential_template, data = data
        )
      ))
    }
  } else {
    return(list(Q = Q))
  }
}

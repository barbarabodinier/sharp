#' Stability selection in Structural Equation Modelling
#'
#' Performs stability selection for Structural Equation Models. The underlying
#' arrow selection algorithm (e.g. regularised Structural Equation Modelling) is
#' run with different combinations of parameters controlling the sparsity (e.g.
#' penalty parameter) and thresholds in selection proportions. These two
#' hyper-parameters are jointly calibrated by maximisation of the stability
#' score.
#'
#' @inheritParams VariableSelection
#' @param xdata matrix with observations as rows and variables as columns.
#' @param adjacency binary adjacency matrix of the Directed Acyclic Graph
#'   (asymmetric matrix A in Reticular Action Model notation).
#' @param residual_covariance binary and symmetric matrix encoding the nonzero
#'   entries in the residual covariance matrix (symmetric matrix S in Reticular
#'   Action Model notation). By default, this is the identity matrix (no
#'   residual covariance).
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#' @param implementation function to use for variable selection. If
#'   \code{implementation=PenalisedSEM}, the algorithm implemented in
#'   \code{\link[regsem]{regsem}} is used.
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
#'   In Structural Equation Modelling, "feature" refers to an arrow in the
#'   corresponding Directed Acyclic Graph.
#'
#'   These parameters can be calibrated by maximisation of a stability score
#'   (see \code{\link{StabilityScore}}) derived from the likelihood under the
#'   assumption of uniform (uninformative) selection:
#'
#'   \eqn{S_{\lambda, \pi} = -log(L_{\lambda, \pi})}
#'
#'   It is strongly recommended to examine the calibration plot carefully to
#'   check that the grids of parameters \code{Lambda} and \code{pi_list} do not
#'   restrict the calibration to a region that would not include the global
#'   maximum (see \code{\link{CalibrationPlot}}). In particular, the grid
#'   \code{Lambda} may need to be extended when the maximum stability is
#'   observed on the left or right edges of the calibration heatmap.
#'
#'   To control the expected number of False Positives (Per Family Error Rate)
#'   in the results, a threshold \code{PFER_thr} can be specified. The
#'   optimisation problem is then constrained to sets of parameters that
#'   generate models with an upper-bound in PFER below \code{PFER_thr} (see
#'   Meinshausen and Bühlmann (2010) and Shah and Samworth (2013)).
#'
#'   Possible resampling procedures include defining (i) \code{K} subsamples of
#'   a proportion \code{tau} of the observations, (ii) \code{K} bootstrap
#'   samples with the full sample size (obtained with replacement), and (iii)
#'   \code{K/2} splits of the data in half for complementary pair stability
#'   selection (see arguments \code{resampling} and \code{cpss}). In
#'   complementary pair stability selection, a feature is considered selected at
#'   a given resampling iteration if it is selected in the two complementary
#'   subsamples.
#'
#'   To ensure reproducibility of the results, the starting number of the random
#'   number generator is set to \code{seed}.
#'
#'   For parallelisation, stability selection with different sets of parameters
#'   can be run on \code{n_cores} cores. This relies on forking with
#'   \code{\link[parallel]{mclapply}} (specific to Unix systems). Alternatively,
#'   the function can be run manually with different \code{seed}s and all other
#'   parameters equal. The results can then be combined using
#'   \code{\link{Combine}}.
#'
#' @return An object of class \code{variable_selection}. A list with: \item{S}{a
#'   matrix of the best stability scores for different parameters controlling
#'   the level of sparsity in the underlying algorithm.} \item{Lambda}{a matrix
#'   of parameters controlling the level of sparsity in the underlying
#'   algorithm.} \item{Q}{a matrix of the average number of selected features by
#'   the underlying algorithm with different parameters controlling the level of
#'   sparsity.} \item{Q_s}{a matrix of the calibrated number of stably selected
#'   features with different parameters controlling the level of sparsity.}
#'   \item{P}{a matrix of calibrated thresholds in selection proportions for
#'   different parameters controlling the level of sparsity in the underlying
#'   algorithm.} \item{PFER}{a matrix of upper-bounds in PFER of calibrated
#'   stability selection models with different parameters controlling the level
#'   of sparsity.} \item{FDP}{a matrix of upper-bounds in FDP of calibrated
#'   stability selection models with different parameters controlling the level
#'   of sparsity.} \item{S_2d}{a matrix of stability scores obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in selection proportions.} \item{PFER_2d}{a matrix of
#'   upper-bounds in FDP obtained with different combinations of parameters.
#'   Columns correspond to different thresholds in selection proportions.}
#'   \item{FDP_2d}{a matrix of upper-bounds in PFER obtained with different
#'   combinations of parameters. Columns correspond to different thresholds in
#'   selection proportions.} \item{selprop}{a matrix of selection proportions.
#'   Columns correspond to predictors from \code{xdata}.} \item{Beta}{an array
#'   of model coefficients. Columns correspond to predictors from \code{xdata}.
#'   Indices along the third dimension correspond to different resampling
#'   iterations. With multivariate outcomes, indices along the fourth dimension
#'   correspond to outcome-specific coefficients.} \item{method}{a list with
#'   \code{type="variable_selection"} and values used for arguments
#'   \code{implementation}, \code{family}, \code{resampling}, \code{cpss} and
#'   \code{PFER_method}.} \item{params}{a list with values used for arguments
#'   \code{K}, \code{pi_list}, \code{tau}, \code{n_cat}, \code{pk}, \code{n}
#'   (number of observations), \code{PFER_thr}, \code{FDP_thr} and \code{seed}.
#'   The datasets \code{xdata} and \code{ydata} are also included if
#'   \code{output_data=TRUE}.} For all matrices and arrays returned, the rows
#'   are ordered in the same way and correspond to parameter values stored in
#'   \code{Lambda}.
#'
#' @family stability selection functions
#' @seealso \code{\link{PenalisedSEM}}, \code{\link{SelectionAlgo}},
#'   \code{\link{Resample}}, \code{\link{StabilityScore}}
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#'   \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{stabilityselectionSS}{sharp}
#'
#'   \insertRef{RegSEM}{sharp}
#'
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = rep(7, 4))
#'
#' # Definition of the model structure
#' layers <- list(
#'   c("var1", "var2", "var3"),
#'   c("var4", "var5"),
#'   c("var6", "var7", "var8")
#' )
#' dag <- LayeredDAG(layers)
#'
#' # Definition of simulated effects
#' theta <- dag
#' theta[2, 4] <- 0
#' theta[3, 7] <- 0
#' theta[4, 7] <- 0
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateStructural(n = 500, theta = theta)
#'
#' # Stability selection
#' stab <- StructuralEquations(
#'   xdata = simul$data,
#'   Lambda = LambdaSequence(lmax = 1, lmin = 0.2, cardinal = 10),
#'   adjacency = dag
#' )
#' CalibrationPlot(stab)
#'
#' # Extracting stable refined DAG
#' LavaanMatrix(SelectedVariables(stab), adjacency = dag)
#' theta # simulated is the same
#'
#' par(oldpar)
#' }
#'
#' @export
StructuralEquations <- function(xdata, adjacency, residual_covariance = NULL,
                                Lambda, pi_list = seq(0.6, 0.9, by = 0.01),
                                K = 100, tau = 0.5, seed = 1, n_cat = 3,
                                implementation = PenalisedSEM,
                                resampling = "subsampling", cpss = FALSE,
                                PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                                Lambda_cardinal = 100,
                                n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  # Object preparation, error and warning messages
  family <- "gaussian"
  ydata <- NULL
  # Object preparation, error and warning messages
  CheckParamRegression(
    Lambda = Lambda, pi_list = pi_list,
    K = K, tau = tau, seed = seed, n_cat = n_cat,
    family = family, implementation = implementation,
    resampling = resampling, PFER_method = PFER_method,
    PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    Lambda_cardinal = Lambda_cardinal,
    verbose = verbose
  )
  CheckDataRegression(
    xdata = xdata, ydata = ydata, family = family, verbose = verbose
  )
  # CheckInputRegression(
  #   xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
  #   K = K, tau = tau, seed = seed, n_cat = n_cat,
  #   family = family, implementation = implementation,
  #   resampling = resampling, PFER_method = PFER_method,
  #   PFER_thr = PFER_thr, FDP_thr = FDP_thr,
  #   Lambda_cardinal = Lambda_cardinal,
  #   verbose = verbose
  # )

  # if (is.null(Lambda)) {
  #   # Defining grid of lambda values (using glmnet implementation)
  #   Lambda <- LambdaGridRegression(
  #     xdata = xdata, ydata = ydata, tau = tau, seed = seed,
  #     family = family,
  #     resampling = resampling,
  #     Lambda_cardinal = Lambda_cardinal, check_input = FALSE, ...
  #   )
  # }

  # Check if parallelisation is possible (forking)
  if (.Platform$OS.type != "unix") {
    if (n_cores > 1) {
      warning("Invalid input for argument 'n_cores'. Parallelisation relies on forking, it is only available on Unix systems.")
    }
    n_cores <- 1
  }

  # Stability selection and score
  mypar <- parallel::mclapply(X = 1:n_cores, FUN = function(k) {
    return(SerialRegression(
      xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
      K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
      family = family, implementation = implementation,
      resampling = resampling, cpss = cpss,
      PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      group_x = NULL, group_penalisation = FALSE,
      output_data = output_data, verbose = verbose,
      adjacency = adjacency, residual_covariance = residual_covariance, ...
    ))
  })

  # Combining the outputs from parallel iterations
  out <- mypar[[1]]
  if (n_cores > 1) {
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[2]]))
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
  }

  # Defining the class
  class(out) <- "variable_selection"

  return(out)
}

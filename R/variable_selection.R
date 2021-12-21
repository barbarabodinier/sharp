#' Stability selection in regression
#'
#' Runs stability selection regression models with different combinations of
#' parameters controlling the sparsity of the underlying selection algorithm
#' (e.g. penalty parameter for regularised models) and thresholds in selection
#' proportions. These two parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features).
#'
#' @param xdata matrix of predictors with observations as rows and variables as
#'   columns.
#' @param ydata optional vector or matrix of outcome(s). If \code{family} is set
#'   to \code{"binomial"} or \code{"multinomial"}, \code{ydata} can be a vector
#'   with character/numeric values, or a factor.
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   With \code{implementation=PenalisedRegression}, \code{Lambda} contains
#'   penalty parameters. If \code{Lambda=NULL},
#'   \code{\link{LambdaGridRegression}} is used to define a relevant grid. If
#'   \code{implementation} is not set to \code{PenalisedRegression},
#'   \code{Lambda} must be provided.
#' @param pi_list vector of thresholds in selection proportions. If
#'   \code{n_cat=3}, these values must be \code{>0.5} and \code{<1}. If
#'   \code{n_cat=2}, these values must be \code{>0} and \code{<1}.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with \code{resampling="subsampling"}.
#' @param seed value of the seed to ensure reproducibility of the results.
#' @param n_cat number of categories used to compute the stability score.
#'   Possible values are 2 or 3.
#' @param family type of regression model. If
#'   \code{implementation=PenalisedRegression}, this argument is defined as in
#'   \code{\link[glmnet]{glmnet}}. Possible values include \code{"gaussian"}
#'   (linear regression), \code{"binomial"} (logistic regression),
#'   \code{"multinomial"} (multinomial regression), and \code{"cox"} (survival
#'   analysis).
#' @param implementation function to use for variable selection. By default,
#'   \code{PenalisedRegression}, based on \code{\link[glmnet]{glmnet}}, is used
#'   for regularised regression. Other possible functions are: \code{SparsePLS},
#'   \code{GroupPLS} and \code{SparseGroupPLS}. Alternatively, a function with
#'   arguments \code{xdata}, \code{ydata}, \code{Lambda}, \code{family} and
#'   \code{...}, and returning a list of two matrices named \code{selected} and
#'   \code{beta_full} of the correct dimensions can be used.
#' @param resampling resampling approach. Possible values are:
#'   \code{"subsampling"} for sampling without replacement of a proportion
#'   \code{tau} of the observations, or \code{"bootstrap"} for sampling with
#'   replacement generating a resampled dataset with as many observations as in
#'   the full sample. Alternatively, this argument can be a function to use for
#'   resampling. This function must use arguments named \code{data} and
#'   \code{tau} and return IDs of observations to be included in the resampled
#'   dataset (see example in \code{\link{Resample}}).
#' @param PFER_method method used to compute the upper-bound of the expected
#'   number of False Positives (or Per Family Error Rate, PFER). With
#'   \code{PFER_method="MB"}, the method proposed by Meinshausen and Bühlmann
#'   (2010) is used. With \code{PFER_method="SS"}, the method proposed by Shah
#'   and Samworth (2013) under the assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error
#'   control. With \code{PFER_thr=Inf} and \code{FDP_thr=Inf}, unconstrained
#'   calibration is used.
#' @param FDP_thr threshold in the expected proportion of falsely selected
#'   features (or False Discovery Proportion, FDP) for constrained calibration
#'   by error control. With \code{PFER_thr=Inf} and \code{FDP_thr=Inf},
#'   unconstrained calibration is used.
#' @param Lambda_cardinal number of values in the grid of parameters controlling
#'   the level of sparsity in the underlying algorithm.
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group. Only used with
#'   \code{implementation=SparseGroupPLS} or \code{implementation=GroupPLS}.
#' @param group_penalisation logical indicating if a group penalisation should
#'   be considered in the stability score. An extra argument \code{group_x}, a
#'   vector encoding the number of variables in each group, must be provided if
#'   \code{group_penalisation=TRUE}. The use of \code{group_penalisation=TRUE}
#'   strictly applies to group (not sparse-group) penalisation.
#' @param n_cores number of cores to use for parallel computing. Only available
#'   on Unix systems.
#' @param output_data logical indicating if the input datasets \code{xdata} and
#'   \code{ydata} should be included in the output.
#' @param verbose logical indicating if a loading bar and messages should be
#'   printed.
#' @param ... additional parameters passed to the functions provided in
#'   \code{implementation} or \code{resampling}.
#'
#' @details To ensure reproducibility of the results, the state of the random
#'   number generator is fixed to \code{seed}. For parallelisation of the code,
#'   stability selection results produced with different \code{seed}s and all
#'   other parameters equal can be combined (more details in
#'   \code{\link{Combine}}).
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
#'   to outcome-specific coefficients.} \item{method}{a list with
#'   \code{type="variable_selection"}, \code{implementation}, \code{family},
#'   \code{resampling} and \code{PFER_method} values used for the run.}
#'   \item{params}{a list of \code{K}, \code{pi_list}, \code{tau}, \code{n_cat},
#'   \code{pk}, \code{n} (number of observations), \code{PFER_thr},
#'   \code{FDP_thr} and \code{seed} values used for the run. The datasets
#'   \code{xdata} and \code{ydata} are also included if
#'   \code{output_data=TRUE}.} For all objects except those stored in
#'   \code{methods} or \code{params}, rows correspond to parameter values stored
#'   in the output \code{Lambda}.
#'
#' @family stability selection functions
#' @seealso \code{\link{Recalibrate}}, \code{\link{ExplanatoryPerformance}},
#'   \code{\link{PlotROC}}, \code{\link{Incremental}},
#'   \code{\link{PlotIncremental}}, \code{\link{Combine}},
#'   \code{\link{LambdaGridRegression}}, \code{\link{Resample}},
#'   \code{\link{StabilityScore}}
#'
#' @references \insertRef{ourstabilityselection}{focus}
#'
#'   \insertRef{stabilityselectionMB}{focus}
#'
#'   \insertRef{stabilityselectionSS}{focus}
#'
#' @examples
#' \dontshow{
#' # Linear regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = 10, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian", K = 5, verbose = FALSE)
#' CalibrationPlot(stab)
#' myselected <- SelectedVariables(stab)
#' perf <- SelectionPerformance(theta = myselected, theta_star = simul$theta[, 1])
#' SelectionProportions(stab)
#' }
#' \dontrun{
#' par(mar = rep(7, 4))
#'
#' # Linear regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' print(stab)
#' CalibrationPlot(stab)
#' summary(stab)
#' SelectedVariables(stab)
#'
#' # Using additional arguments from glmnet (e.g. penalty.factor)
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata, family = "gaussian",
#'   penalty.factor = c(rep(1, 45), rep(0, 5))
#' )
#' summary(stab)
#'
#' # Regression with multivariate outcomes
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = c(20, 30), family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "mgaussian")
#' summary(stab)
#'
#' # Logistic regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "binomial")
#' summary(stab)
#'
#' # Multinomial regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 15, family = "multinomial")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "multinomial"
#' )
#' summary(stab)
#'
#' # Sparse PCA (1 component, see BiSelection for more components)
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4))
#' stab <- VariableSelection(
#'   xdata = simul$data,
#'   Lambda = 1:(ncol(simul$data) - 1),
#'   implementation = SparsePCA
#' )
#' CalibrationPlot(stab, xlab = "")
#' summary(stab)
#'
#' # Sparse PLS (1 outcome, 1 component, see BiSelection for more options)
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = 1:(ncol(simul$xdata) - 1),
#'   implementation = SparsePLS, family = "gaussian"
#' )
#' CalibrationPlot(stab, xlab = "")
#' SelectedVariables(stab)
#'
#' # Group PLS (1 outcome, 1 component, see BiSelection for more options)
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = 1:5,
#'   group_x = c(5, 5, 10, 20, 10),
#'   group_penalisation = TRUE,
#'   implementation = GroupPLS, family = "gaussian"
#' )
#' CalibrationPlot(stab, xlab = "")
#' SelectedVariables(stab)
#'
#' # Sparse PLS-DA (1 outcome, 1 component, see BiSelection for more options)
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   Lambda = 1:(ncol(simul$xdata) - 1),
#'   implementation = SparsePLS,
#'   family = "binomial"
#' )
#' CalibrationPlot(stab, xlab = "")
#' summary(stab)
#'
#' # Example with more hyper-parameters: elastic net
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' TuneElasticNet <- function(xdata, ydata, family, alpha) {
#'   stab <- VariableSelection(
#'     xdata = xdata, ydata = ydata,
#'     family = family, alpha = alpha, verbose = FALSE
#'   )
#'   return(max(stab$S, na.rm = TRUE))
#' }
#' myopt <- optimise(TuneElasticNet,
#'   lower = 0.1, upper = 1, maximum = TRUE,
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian"
#' )
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", alpha = myopt$maximum
#' )
#' summary(stab)
#' enet <- SelectedVariables(stab)
#'
#' # Comparison with LASSO
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' summary(stab)
#' lasso <- SelectedVariables(stab)
#' table(lasso, enet)
#'
#' # Example using an external function: group-LASSO with gglasso
#' if (requireNamespace("gglasso", quietly = TRUE)) {
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#'   ManualGridGroupLasso <- function(xdata, ydata, family, group_x, ...) {
#'     # Defining the grouping
#'     group <- do.call(c, lapply(1:length(group_x), FUN = function(i) {
#'       rep(i, group_x[i])
#'     }))
#'
#'     if (family == "binomial") {
#'       ytmp <- ydata
#'       ytmp[ytmp == min(ytmp)] <- -1
#'       ytmp[ytmp == max(ytmp)] <- 1
#'       return(gglasso::gglasso(xdata, ytmp, loss = "logit", group = group, ...))
#'     } else {
#'       return(gglasso::gglasso(xdata, ydata, lambda = lambda, loss = "ls", group = group, ...))
#'     }
#'   }
#'   Lambda <- LambdaGridRegression(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     family = "binomial", Lambda_cardinal = 20,
#'     implementation = ManualGridGroupLasso,
#'     group_x = rep(5, 4)
#'   )
#'   GroupLasso <- function(xdata, ydata, Lambda, family, group_x, ...) {
#'     # Defining the grouping
#'     group <- do.call(c, lapply(1:length(group_x), FUN = function(i) {
#'       rep(i, group_x[i])
#'     }))
#'
#'     # Running the regression
#'     if (family == "binomial") {
#'       ytmp <- ydata
#'       ytmp[ytmp == min(ytmp)] <- -1
#'       ytmp[ytmp == max(ytmp)] <- 1
#'       mymodel <- gglasso::gglasso(xdata, ytmp, lambda = Lambda, loss = "logit", group = group, ...)
#'     }
#'     if (family == "gaussian") {
#'       mymodel <- gglasso::gglasso(xdata, ydata, lambda = Lambda, loss = "ls", group = group, ...)
#'     }
#'     # Extracting and formatting the beta coefficients
#'     beta_full <- t(as.matrix(mymodel$beta))
#'     beta_full <- beta_full[, colnames(xdata)]
#'
#'     selected <- ifelse(beta_full != 0, yes = 1, no = 0)
#'
#'     return(list(selected = selected, beta_full = beta_full))
#'   }
#'   stab <- VariableSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     implementation = GroupLasso, family = "binomial", Lambda = Lambda,
#'     group_x = rep(5, 4),
#'     group_penalisation = TRUE
#'   )
#'   summary(stab)
#' }
#' }
#' @export
VariableSelection <- function(xdata, ydata = NULL, Lambda = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                              K = 100, tau = 0.5, seed = 1, n_cat = 3,
                              family = "gaussian", implementation = PenalisedRegression,
                              resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                              Lambda_cardinal = 100, group_x = NULL, group_penalisation = FALSE,
                              n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  # Defining Lambda if used with sparse PCA or PLS
  if (is.null(Lambda)) {
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", "GroupPLS")) {
      Lambda <- seq(1, length(group_x) - 1)
    }
    if (as.character(substitute(implementation)) %in% c("SparsePLS", "SparsePCA")) {
      Lambda <- seq(1, ncol(xdata) - 1)
    }
  }

  # Object preparation, error and warning messages
  CheckInputRegression(
    xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
    K = K, tau = tau, seed = seed, n_cat = n_cat,
    family = family, implementation = implementation,
    resampling = resampling, PFER_method = PFER_method,
    PFER_thr = PFER_thr, FDP_thr = FDP_thr,
    Lambda_cardinal = Lambda_cardinal,
    verbose = verbose
  )

  # Checking that group_x is provided for group penalisation
  if (group_penalisation) {
    if (is.null(group_x)) {
      stop("Please provide argument 'group_x' for group penalisation. Argument 'group_x' should be a vector with the number of variables in each group.")
    }
  }

  if (is.null(Lambda)) {
    # Defining grid of lambda values (using glmnet implementation)
    Lambda <- LambdaGridRegression(
      xdata = xdata, ydata = ydata, tau = tau, seed = seed,
      family = family,
      resampling = resampling,
      Lambda_cardinal = Lambda_cardinal, check_input = FALSE, ...
    )
  }

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
      family = family, implementation = implementation, resampling = resampling,
      PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      group_x = group_x, group_penalisation = group_penalisation,
      output_data = output_data, verbose = verbose, ...
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


#' Stability selection in regression (internal)
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
#' @inheritParams VariableSelection
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
SerialRegression <- function(xdata, ydata = NULL, Lambda, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, tau = 0.5, seed = 1, n_cat = 3,
                             family = "gaussian", implementation = PenalisedRegression,
                             resampling = "subsampling", PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                             group_x = NULL, group_penalisation = FALSE,
                             output_data = FALSE, verbose = TRUE, ...) {
  # Defining K if using complementary pairs (SS)
  if (PFER_method == "SS") {
    K <- ceiling(K / 2) * 2
    tau <- 0.5
  }

  # Initialising objects to be filled
  N <- N_block <- ncol(xdata)

  # Initialising the arrays
  s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
  Xsub <- xdata[s, ]
  Ysub <- ydata[s, ]
  mybeta <- SelectionAlgo(
    xdata = Xsub, ydata = Ysub,
    Lambda = Lambda[, 1], group_x = group_x,
    family = family, implementation = implementation, ...
  )
  Beta <- array(0, dim = c(nrow(mybeta$selected), ncol(mybeta$selected), K))
  rownames(Beta) <- rownames(mybeta$selected)
  colnames(Beta) <- colnames(mybeta$selected)
  if (length(dim(mybeta$beta_full)) == 2) {
    Beta_full <- array(0,
      dim = c(nrow(Lambda), dim(mybeta$beta_full)[2], K),
      dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL)
    )
  } else {
    if (length(dim(mybeta$beta_full)) == 3) {
      Beta_full <- array(0,
        dim = c(nrow(Lambda), dim(mybeta$beta_full)[2], K, dim(mybeta$beta_full)[3]),
        dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL, dimnames(mybeta$beta_full)[[3]])
      )
    } else {
      stop(paste0("Invalid output from the variable selection function: ", implementation, "(). The output 'beta_full' must be an array with 2 or 3 dimensions."))
    }
  }

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Computation of the selection proportions over Lambda
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  if (PFER_method == "MB") {
    for (k in 1:K) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
      Xsub <- xdata[s, ]
      Ysub <- ydata[s, ]
      mybeta <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Resampling if model failed to converge
      while (is.infinite(mybeta$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
        Xsub <- xdata[s, ]
        Ysub <- ydata[s, ]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )
      }

      # Storing (one set of) beta coefficients, used to define set of selected variables
      Beta[rownames(mybeta$selected), colnames(mybeta$selected), k] <- mybeta$selected

      # Storing all beta coefficients
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta$beta_full), colnames(mybeta$beta_full), k] <- mybeta$beta_full
      } else {
        Beta_full[rownames(mybeta$beta_full), colnames(mybeta$beta_full), k, ] <- mybeta$beta_full
      }

      if (verbose) {
        utils::setTxtProgressBar(pb, k / K)
      }
    }

    # Computing the selection proportions
    bigstab <- matrix(NA, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab) <- colnames(Beta)
    rownames(bigstab) <- rownames(Beta)
    for (i in 1:nrow(Beta)) {
      for (j in 1:ncol(Beta)) {
        bigstab[i, j] <- sum(Beta[i, j, ] != 0) / K
      }
    }
  }

  if (PFER_method == "SS") {
    for (k in 1:ceiling(K / 2)) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)

      # First subset
      Xsub <- xdata[s, ]
      Ysub <- ydata[s, ]
      mybeta1 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Complementary subset
      Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], ]
      Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], ]
      mybeta2 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Resampling if model failed to converge
      while (is.infinite(mybeta1$selected[1]) | is.infinite(mybeta2$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)

        # First subset
        Xsub <- xdata[s, ]
        Ysub <- ydata[s, ]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )

        # Complementary subset
        Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], ]
        Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], ]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )
      }

      # Storing beta coefficients from first set
      Beta[rownames(mybeta1$selected), colnames(mybeta1$selected), k] <- mybeta1$selected

      # Storing all beta coefficients from first set
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta1$beta_full), colnames(mybeta1$beta_full), k] <- mybeta1$beta_full
      } else {
        Beta_full[rownames(mybeta1$beta_full), colnames(mybeta1$beta_full), k, ] <- mybeta1$beta_full
      }

      # Storing beta coefficients from complementary set
      Beta[rownames(mybeta2$selected), colnames(mybeta2$selected), ceiling(K / 2) + k] <- mybeta2$selected

      # Storing all beta coefficients from complementary set
      if (length(dim(Beta_full)) == 3) {
        Beta_full[rownames(mybeta2$beta_full), colnames(mybeta2$beta_full), ceiling(K / 2) + k] <- mybeta2$beta_full
      } else {
        Beta_full[rownames(mybeta2$beta_full), colnames(mybeta2$beta_full), ceiling(K / 2) + k, ] <- mybeta2$beta_full
      }

      if (verbose) {
        utils::setTxtProgressBar(pb, 2 * k / K)
      }
    }

    # Computing the simultaneous selection proportions
    bigstab <- matrix(0, nrow = nrow(Beta), ncol = ncol(Beta))
    colnames(bigstab) <- colnames(Beta)
    rownames(bigstab) <- rownames(Beta)
    for (k in 1:ceiling(K / 2)) {
      A1 <- ifelse(Beta[, , k] != 0, yes = 1, no = 0)
      A2 <- ifelse(Beta[, , ceiling(K / 2) + k] != 0, yes = 1, no = 0)
      A <- A1 + A2
      A <- ifelse(A == 2, yes = 1, no = 0)
      bigstab <- bigstab + A
    }
    bigstab <- bigstab / ceiling(K / 2)
  }

  if (verbose) {
    cat("\n")
  }

  # Computation of the stability score over Lambda and pi_list
  if (group_penalisation) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = NULL, graph = FALSE, group = group_x,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr
    )
  } else {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = NULL, graph = FALSE,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr
    )
  }
  if (verbose) {
    utils::setTxtProgressBar(pb, 1)
    cat("\n")
  }
  Beta <- Beta_full

  # Preparing outputs
  myimplementation <- as.character(substitute(implementation, env = parent.frame(n = 2)))
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  out <- list(
    S = metrics$S, Lambda = Lambda,
    Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
    PFER = metrics$PFER, FDP = metrics$FDP,
    S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
    selprop = bigstab, Beta = Beta,
    methods = list(
      type = "variable_selection", implementation = myimplementation, family = family,
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
    out$params <- c(out$params, list(xdata = xdata, ydata = ydata))
  }

  return(out)
}

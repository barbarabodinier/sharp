#' Stability selection in regression
#'
#' Performs stability selection for regression models. The underlying variable
#' selection algorithm (e.g. LASSO regression) is run with different
#' combinations of parameters controlling the sparsity (e.g. penalty parameter)
#' and thresholds in selection proportions. These two hyper-parameters are
#' jointly calibrated by maximisation of the stability score.
#'
#' @param xdata matrix of predictors with observations as rows and variables as
#'   columns.
#' @param ydata optional vector or matrix of outcome(s). If \code{family} is set
#'   to \code{"binomial"} or \code{"multinomial"}, \code{ydata} can be a vector
#'   with character/numeric values or a factor.
#' @param Lambda matrix of parameters controlling the level of sparsity in the
#'   underlying feature selection algorithm specified in \code{implementation}.
#'   If \code{Lambda=NULL} and \code{implementation=PenalisedRegression},
#'   \code{\link{LambdaGridRegression}} is used to define a relevant grid.
#' @param pi_list vector of thresholds in selection proportions. If
#'   \code{n_cat=NULL} or \code{n_cat=2}, these values must be \code{>0} and
#'   \code{<1}. If \code{n_cat=3}, these values must be \code{>0.5} and
#'   \code{<1}.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used if \code{resampling="subsampling"} and
#'   \code{cpss=FALSE}.
#' @param seed value of the seed to initialise the random number generator and
#'   ensure reproducibility of the results (see \code{\link[base]{set.seed}}).
#' @param n_cat computation options for the stability score. Default is
#'   \code{NULL} to use the score based on a z test. Other possible values are 2
#'   or 3 to use the score based on the negative log-likelihood.
#' @param family type of regression model. This argument is defined as in
#'   \code{\link[glmnet]{glmnet}}. Possible values include \code{"gaussian"}
#'   (linear regression), \code{"binomial"} (logistic regression),
#'   \code{"multinomial"} (multinomial regression), and \code{"cox"} (survival
#'   analysis).
#' @param implementation function to use for variable selection. Possible
#'   functions are: \code{PenalisedRegression}, \code{SparsePLS},
#'   \code{GroupPLS} and \code{SparseGroupPLS}. Alternatively, a user-defined
#'   function can be provided.
#' @param resampling resampling approach. Possible values are:
#'   \code{"subsampling"} for sampling without replacement of a proportion
#'   \code{tau} of the observations, or \code{"bootstrap"} for sampling with
#'   replacement generating a resampled dataset with as many observations as in
#'   the full sample. Alternatively, this argument can be a function to use for
#'   resampling. This function must use arguments named \code{data} and
#'   \code{tau} and return the IDs of observations to be included in the
#'   resampled dataset.
#' @param cpss logical indicating if complementary pair stability selection
#'   should be done. For this, the algorithm is applied on two non-overlapping
#'   subsets of half of the observations. A feature is considered as selected if
#'   it is selected for both subsamples. With this method, the data is split
#'   \code{K/2} times (\code{K} models are fitted). Only used if
#'   \code{PFER_method="MB"}.
#' @param PFER_method method used to compute the upper-bound of the expected
#'   number of False Positives (or Per Family Error Rate, PFER). If
#'   \code{PFER_method="MB"}, the method proposed by Meinshausen and Bühlmann
#'   (2010) is used. If \code{PFER_method="SS"}, the method proposed by Shah and
#'   Samworth (2013) under the assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error
#'   control. If \code{PFER_thr=Inf} and \code{FDP_thr=Inf}, unconstrained
#'   calibration is used (the default).
#' @param FDP_thr threshold in the expected proportion of falsely selected
#'   features (or False Discovery Proportion) for constrained calibration by
#'   error control. If \code{PFER_thr=Inf} and \code{FDP_thr=Inf}, unconstrained
#'   calibration is used (the default).
#' @param Lambda_cardinal number of values in the grid of parameters controlling
#'   the level of sparsity in the underlying algorithm. Only used if
#'   \code{Lambda=NULL}.
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group. Only used for
#'   models with group penalisation (e.g. \code{implementation=GroupPLS} or
#'   \code{implementation=SparseGroupPLS}).
#' @param group_penalisation logical indicating if a group penalisation should
#'   be considered in the stability score. The use of
#'   \code{group_penalisation=TRUE} strictly applies to group (not sparse-group)
#'   penalisation.
#' @param optimisation character string indicating the type of optimisation
#'   method. With \code{optimisation="grid_search"} (the default), all values in
#'   \code{Lambda} are visited. Alternatively, optimisation algorithms
#'   implemented in \code{\link[nloptr]{nloptr}} can be used with
#'   \code{optimisation="nloptr"}. By default, we use
#'   \code{"algorithm"="NLOPT_GN_DIRECT_L"}, \code{"xtol_abs"=0.1},
#'   \code{"ftol_abs"=0.1} and \code{"maxeval"=Lambda_cardinal}. These values
#'   can be changed by providing the argument \code{opts} (see
#'   \code{\link[nloptr]{nloptr}}). For stability selection using penalised
#'   regression, \code{optimisation="grid_search"} may be faster as it allows
#'   for warm start.
#' @param n_cores number of cores to use for parallel computing (see argument
#'   \code{workers} in \code{\link[future]{multisession}}). Using
#'   \code{n_cores>1} is only supported with \code{optimisation="grid_search"}.
#' @param output_data logical indicating if the input datasets \code{xdata} and
#'   \code{ydata} should be included in the output.
#' @param verbose logical indicating if a loading bar and messages should be
#'   printed.
#' @param beep sound indicating the end of the run. Possible values are:
#'   \code{NULL} (no sound) or an integer between 1 and 11 (see argument
#'   \code{sound} in \code{\link[beepr]{beep}}).
#' @param ... additional parameters passed to the functions provided in
#'   \code{implementation} or \code{resampling}.
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
#'   If argument \code{group_penalisation=FALSE}, "feature" refers to variable
#'   (variable selection model). If argument \code{group_penalisation=TRUE},
#'   "feature" refers to group (group selection model). In this case, groups
#'   need to be defined \emph{a priori} and specified in argument
#'   \code{group_x}.
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
#'   a proportion \code{tau} of the observations, (ii) \code{K} bootstrap
#'   samples with the full sample size (obtained with replacement), and (iii)
#'   \code{K/2} splits of the data in half for complementary pair stability
#'   selection (see arguments \code{resampling} and \code{cpss}). In
#'   complementary pair stability selection, a feature is considered selected at
#'   a given resampling iteration if it is selected in the two complementary
#'   subsamples.
#'
#'   For categorical or time to event outcomes (argument \code{family} is
#'   \code{"binomial"}, \code{"multinomial"} or \code{"cox"}), the proportions
#'   of observations from each category in all subsamples or bootstrap samples
#'   are the same as in the full sample.
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
#' @family stability functions
#' @seealso \code{\link{PenalisedRegression}}, \code{\link{SelectionAlgo}},
#'   \code{\link{LambdaGridRegression}}, \code{\link{Resample}},
#'   \code{\link{StabilityScore}} \code{\link{Refit}},
#'   \code{\link{ExplanatoryPerformance}}, \code{\link{Incremental}},
#'
#' @references \insertRef{JStatSoft}{sharp}
#'
#'   \insertRef{ourstabilityselection}{sharp}
#'
#'   \insertRef{stabilityselectionSS}{sharp}
#'
#'   \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{lasso}{sharp}
#'
#' @import fake
#'
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = rep(7, 4))
#'
#' # Linear regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian"
#' )
#'
#' # Calibration plot
#' CalibrationPlot(stab)
#'
#' # Extracting the results
#' summary(stab)
#' Stable(stab)
#' SelectionProportions(stab)
#' plot(stab)
#'
#' # Using randomised LASSO
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", penalisation = "randomised"
#' )
#' plot(stab)
#'
#' # Using adaptive LASSO
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", penalisation = "adaptive"
#' )
#' plot(stab)
#'
#' # Using additional arguments from glmnet (e.g. penalty.factor)
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata, family = "gaussian",
#'   penalty.factor = c(rep(1, 45), rep(0, 5))
#' )
#' head(coef(stab))
#'
#' # Using CART
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   stab <- VariableSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     implementation = CART,
#'     family = "gaussian",
#'   )
#'   plot(stab)
#' }
#'
#' # Regression with multivariate outcomes
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 20, q = 3, family = "gaussian")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "mgaussian"
#' )
#' summary(stab)
#'
#' # Logistic regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 10, family = "binomial", ev_xy = 0.8)
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "binomial"
#' )
#' summary(stab)
#'
#' # Sparse PCA (1 component, see BiSelection for more components)
#' if (requireNamespace("elasticnet", quietly = TRUE)) {
#'   set.seed(1)
#'   simul <- SimulateComponents(pk = c(5, 3, 4))
#'   stab <- VariableSelection(
#'     xdata = simul$data,
#'     Lambda = seq_len(ncol(simul$data) - 1),
#'     implementation = SparsePCA
#'   )
#'   CalibrationPlot(stab, xlab = "")
#'   summary(stab)
#' }
#'
#' # Sparse PLS (1 outcome, 1 component, see BiSelection for more options)
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'   stab <- VariableSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     Lambda = seq_len(ncol(simul$xdata) - 1),
#'     implementation = SparsePLS, family = "gaussian"
#'   )
#'   CalibrationPlot(stab, xlab = "")
#'   SelectedVariables(stab)
#' }
#'
#' # Group PLS (1 outcome, 1 component, see BiSelection for more options)
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   stab <- VariableSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     Lambda = seq_len(5),
#'     group_x = c(5, 5, 10, 20, 10),
#'     group_penalisation = TRUE,
#'     implementation = GroupPLS, family = "gaussian"
#'   )
#'   CalibrationPlot(stab, xlab = "")
#'   SelectedVariables(stab)
#' }
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
#'     group <- do.call(c, lapply(seq_len(length(group_x)), FUN = function(i) {
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
#'     group <- do.call(c, lapply(seq_len(length(group_x)), FUN = function(i) {
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
#'
#' par(oldpar)
#' }
#'
#' @export
VariableSelection <- function(xdata, ydata = NULL, Lambda = NULL, pi_list = seq(0.01, 0.99, by = 0.01),
                              K = 100, tau = 0.5, seed = 1, n_cat = NULL,
                              family = "gaussian", implementation = PenalisedRegression,
                              resampling = "subsampling", cpss = FALSE,
                              PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                              Lambda_cardinal = 100, group_x = NULL, group_penalisation = FALSE,
                              optimisation = c("grid_search", "nloptr"),
                              n_cores = 1, output_data = FALSE, verbose = TRUE, beep = NULL, ...) {
  if (is.null(Lambda)) {
    # Defining Lambda if used with sparse PLS
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", "GroupPLS")) {
      Lambda <- seq(1, length(group_x) - 1)
    }

    # Defining Lambda if used with sparse PCA
    if (as.character(substitute(implementation)) %in% c("SparsePLS", "SparsePCA")) {
      Lambda <- seq(1, ncol(xdata) - 1)
    }

    # Defining Lambda if used with CART
    if (as.character(substitute(implementation)) %in% c("CART")) {
      Lambda <- seq(1, min(nrow(xdata) / 2, 100))
    }
  }

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

  # Defining the type of optimisation
  optimisation <- match.arg(optimisation)

  # Storing extra arguments
  extra_args <- list(...)

  # Stability selection and score
  if (n_cores > 1) {
    if (optimisation != "grid_search") {
      message("Using grid search to allow for parallelisation.")
    }
    future::plan(future::multisession, workers = n_cores)
    mypar <- future.apply::future_lapply(X = seq_len(n_cores), future.seed = TRUE, FUN = function(k) {
      return(SerialRegression(
        xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
        K = ceiling(K / n_cores), tau = tau, seed = as.numeric(paste0(seed, k)), n_cat = n_cat,
        family = family, implementation = implementation,
        resampling = resampling, cpss = cpss,
        PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
        group_x = group_x, group_penalisation = group_penalisation,
        output_data = output_data, verbose = FALSE, ...
      ))
    })
    future::plan(future::sequential)

    # Combining the outputs from parallel iterations
    out <- mypar[[1]]
    if (n_cores > 1) {
      for (i in 2:length(mypar)) {
        out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[i]]))
      }
    }
  } else {
    if (optimisation == "grid_search") {
      out <- SerialRegression(
        xdata = xdata, ydata = ydata, Lambda = Lambda, pi_list = pi_list,
        K = K, tau = tau, seed = seed, n_cat = n_cat,
        family = family, implementation = implementation,
        resampling = resampling, cpss = cpss,
        PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
        group_x = group_x, group_penalisation = group_penalisation,
        output_data = output_data, verbose = verbose, ...
      )
    } else {
      # Creating the function to be minimised
      eval_f <- function(x, env) {
        # Running with a given lambda
        out_nloptr <- SerialRegression(
          xdata = xdata, ydata = ydata, Lambda = rbind(x), pi_list = pi_list,
          K = K, tau = tau, seed = seed, n_cat = n_cat,
          family = family, implementation = implementation,
          resampling = resampling, cpss = cpss,
          PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr,
          group_x = group_x, group_penalisation = group_penalisation,
          output_data = output_data, verbose = FALSE, ...
        )
        if (any(!is.na(out_nloptr$S))) {
          score <- max(out_nloptr$S, na.rm = TRUE)
        } else {
          score <- -Inf
        }

        # Storing the visited values
        out <- get("out", envir = env)
        out <- Concatenate(out_nloptr, out)
        assign("out", out, envir = env)

        return(-score)
      }

      # Defining the nloptr options
      opts <- list(
        "algorithm" = "NLOPT_GN_DIRECT_L",
        "xtol_abs" = 0.1,
        "ftol_abs" = 0.1,
        "print_level" = 0,
        "maxeval" = Lambda_cardinal
      )
      if ("opts" %in% names(extra_args)) {
        for (opts_id in 1:length(extra_args[["opts"]])) {
          opts[[names(extra_args[["opts"]])[opts_id]]] <- extra_args[["opts"]][[opts_id]]
        }
      }

      # Initialising the values to store
      nloptr_env <- new.env(parent = emptyenv())
      assign("out", NULL, envir = nloptr_env)
      nloptr_results <- nloptr::nloptr(
        x0 = apply(Lambda, 2, max, na.rm = TRUE),
        eval_f = eval_f,
        opts = opts,
        lb = apply(Lambda, 2, min, na.rm = TRUE),
        ub = apply(Lambda, 2, max, na.rm = TRUE),
        env = nloptr_env
      )
      out <- get("out", envir = nloptr_env)
      out <- Concatenate(out, order_output = TRUE)
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

  # Making beep
  if (!is.null(beep)) {
    beepr::beep(sound = beep)
  }

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
#'   number of selected features by the underlying algorithm with different
#'   parameters controlling the level of sparsity.} \item{Q_s}{a matrix of the
#'   calibrated number of stably selected features with different parameters
#'   controlling the level of sparsity.} \item{P}{a matrix of calibrated
#'   thresholds in selection proportions for different parameters controlling
#'   the level of sparsity in the underlying algorithm.} \item{PFER}{a matrix of
#'   upper-bounds in PFER of calibrated stability selection models with
#'   different parameters controlling the level of sparsity.} \item{FDP}{a
#'   matrix of upper-bounds in FDP of calibrated stability selection models with
#'   different parameters controlling the level of sparsity.} \item{S_2d}{a
#'   matrix of stability scores obtained with different combinations of
#'   parameters. Columns correspond to different thresholds in selection
#'   proportions.} \item{PFER_2d}{a matrix of upper-bounds in FDP obtained with
#'   different combinations of parameters. Columns correspond to different
#'   thresholds in selection proportions.} \item{FDP_2d}{a matrix of
#'   upper-bounds in PFER obtained with different combinations of parameters.
#'   Columns correspond to different thresholds in selection proportions.}
#'   \item{selprop}{a matrix of selection proportions. Columns correspond to
#'   predictors from \code{xdata}.} \item{Beta}{an array of model coefficients.
#'   Columns correspond to predictors from \code{xdata}. Indices along the third
#'   dimension correspond to different resampling iterations. With multivariate
#'   outcomes, indices along the fourth dimension correspond to outcome-specific
#'   coefficients.} \item{method}{a list with \code{type="variable_selection"}
#'   and values used for arguments \code{implementation}, \code{family},
#'   \code{resampling}, \code{cpss} and \code{PFER_method}.} \item{params}{a
#'   list with values used for arguments \code{K}, \code{pi_list}, \code{tau},
#'   \code{n_cat}, \code{pk}, \code{n} (number of observations),
#'   \code{PFER_thr}, \code{FDP_thr} and \code{seed}. The datasets \code{xdata}
#'   and \code{ydata} are also included if \code{output_data=TRUE}.} For all
#'   matrices and arrays returned, the rows are ordered in the same way and
#'   correspond to parameter values stored in \code{Lambda}.
#'
#' @keywords internal
SerialRegression <- function(xdata, ydata = NULL, Lambda, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, tau = 0.5, seed = 1, n_cat = 3,
                             family = "gaussian", implementation = PenalisedRegression,
                             resampling = "subsampling", cpss = FALSE,
                             PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf,
                             group_x = NULL, group_penalisation = FALSE,
                             output_data = FALSE, verbose = TRUE, ...) {
  # Using complementary pairs for SS
  if (PFER_method == "SS") {
    cpss <- TRUE
  }

  # Defining K if using complementary pairs
  if (cpss) {
    K <- ceiling(K / 2) * 2
    tau <- 0.5
  }

  # Initialising objects to be filled
  N <- N_block <- ncol(xdata)

  # Initialising the arrays
  s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
  Xsub <- xdata[s, , drop = FALSE]
  Ysub <- ydata[s, , drop = FALSE]
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
      dim = c(dim(mybeta$beta_full)[1], dim(mybeta$beta_full)[2], K),
      dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL)
    )
  } else {
    if (length(dim(mybeta$beta_full)) == 3) {
      Beta_full <- array(0,
        dim = c(dim(mybeta$beta_full)[1], dim(mybeta$beta_full)[2], K, dim(mybeta$beta_full)[3]),
        dimnames = list(rownames(mybeta$beta_full), dimnames(mybeta$beta_full)[[2]], NULL, dimnames(mybeta$beta_full)[[3]])
      )
    } else {
      stop(paste0("Invalid output from the variable selection function provided in 'implementation'. The output 'beta_full' must be an array with 2 or 3 dimensions."))
    }
  }

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Computation of the selection proportions over Lambda
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  if (!cpss) {
    for (k in seq_len(K)) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
      Xsub <- xdata[s, , drop = FALSE]
      Ysub <- ydata[s, , drop = FALSE]
      mybeta <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Resampling if model failed to converge
      while (is.infinite(mybeta$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)
        Xsub <- xdata[s, , drop = FALSE]
        Ysub <- ydata[s, , drop = FALSE]
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
    for (i in seq_len(nrow(Beta))) {
      for (j in seq_len(ncol(Beta))) {
        bigstab[i, j] <- sum(Beta[i, j, ] != 0, na.rm = TRUE) / sum(!is.na(Beta[i, j, ]))
      }
    }
  } else {
    for (k in seq_len(ceiling(K / 2))) {
      s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)

      # First subset
      Xsub <- xdata[s, , drop = FALSE]
      Ysub <- ydata[s, , drop = FALSE]
      mybeta1 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Complementary subset
      Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
      Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
      mybeta2 <- SelectionAlgo(
        xdata = Xsub, ydata = Ysub,
        Lambda = Lambda[, 1], group_x = group_x,
        family = family, implementation = implementation, ...
      )

      # Resampling if model failed to converge
      while (is.infinite(mybeta1$selected[1]) | is.infinite(mybeta2$selected[1])) {
        s <- Resample(data = ydata, family = family, tau = tau, resampling = resampling, ...)

        # First subset
        Xsub <- xdata[s, , drop = FALSE]
        Ysub <- ydata[s, , drop = FALSE]
        mybeta <- SelectionAlgo(
          xdata = Xsub, ydata = Ysub,
          Lambda = Lambda[, 1], group_x = group_x,
          family = family, implementation = implementation, ...
        )

        # Complementary subset
        Xsub <- xdata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
        Ysub <- ydata[seq(1, nrow(xdata))[!seq(1, nrow(xdata)) %in% s], , drop = FALSE]
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
    n_valid <- rep(0, nrow(Beta))
    for (k in seq_len(ceiling(K / 2))) {
      A1 <- ifelse(Beta[, , k] != 0, yes = 1, no = 0)
      A2 <- ifelse(Beta[, , ceiling(K / 2) + k] != 0, yes = 1, no = 0)
      A <- A1 + A2
      A <- ifelse(A == 2, yes = 1, no = 0)
      tmp_missing <- apply(A, 1, FUN = function(x) {
        any(is.na(x))
      })
      A[which(tmp_missing), ] <- 0
      n_valid <- n_valid + ifelse(!tmp_missing, yes = 1, no = 0)
      bigstab <- bigstab + A
    }
    bigstab <- bigstab / n_valid
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
      resampling = myresampling, cpss = cpss, PFER_method = PFER_method
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

  # Defining the class
  class(out) <- "variable_selection"

  return(out)
}

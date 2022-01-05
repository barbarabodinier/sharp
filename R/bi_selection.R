#' Stability selection of both predictors and outcomes
#'
#' Runs stability selection regression models with different combinations of
#' parameters controlling the sparsity in PLS/PCA models and thresholds in
#' selection proportions. These parameters are jointly calibrated by maximising
#' the stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features).
#'
#' @inheritParams VariableSelection
#' @param family type of PLS model. This parameter must be set to
#'   \code{family="gaussian"} for continuous outcomes, or to
#'   \code{family="binomial"} for categorical outcomes. Only used if
#'   \code{ydata} is provided.
#' @param group_y optional vector encoding the grouping structure among
#'   outcomes. This argument indicates the number of variables in each group.
#'   Only used with \code{implementation=SparseGroupPLS}.
#' @param LambdaX matrix of parameters controlling the number of selected
#'   variables (sparse PLS) or groups (sparse group PLS) in X.
#' @param LambdaY matrix of parameters controlling the number of selected
#'   variables (sparse PLS) or groups (sparse group PLS) in Y. Only used with
#'   \code{family="gaussian"}.
#' @param AlphaX matrix of parameters controlling the level of sparsity within
#'   groups (sparse group PLS) in X. Only used with
#'   \code{implementation=SparseGroupPLS}.
#' @param AlphaY matrix of parameters controlling the level of sparsity within
#'   groups (sparse group PLS) in X. Only used with
#'   \code{implementation=SparseGroupPLS} and \code{family="gaussian"}.
#' @param ncomp number of components.
#' @param scale logical indicating if the data should be scaled (i.e.
#'   transformed so that all variables have a standard deviation of one).
#'
#' @return A list with: \item{summary}{a matrix of the best stability scores and
#'   corresponding parameters controlling the level of sparsity in the
#'   underlying algorithm for different numbers of components. Possible columns
#'   include: \code{comp} (component ID), \code{nx} (number of predictors to
#'   include, parameter of the underlying algorithm), \code{alphax} (sparsity
#'   within the predictor groups, parameter of the underlying algorithm),
#'   \code{pix} (threshold in selection proportion for predictors), \code{ny}
#'   (number of outcomes to include, parameter of the underlying algorithm),
#'   \code{alphay} (sparsity within the outcome groups, parameter of the
#'   underlying algorithm), \code{piy} (threshold in selection proportion for
#'   outcomes), \code{S} (stability score). Columns that are not relevant to the
#'   model are not reported. For example, \code{alpha_x} and \code{alpha_y} are
#'   not used in sparse PLS models.} \item{summary_full}{a matrix of the
#'   stability scores for the visited combinations of parameters. Values are
#'   reported for the calibrated threshold(s) in selection proportions.}
#'   \item{selectedX}{a binary matrix encoding stably selected predictors.}
#'   \item{selpropX}{a matrix with the calibrated selection proportions of
#'   predictors.} \item{selectedY}{a binary matrix encoding stably selected
#'   outcomes.} \item{selpropY}{a matrix with the calibrated selection
#'   proportions of outcomes.} \item{selectedX_full}{a binary matrix encoding
#'   stably selected predictors.} \item{selpropX_full}{a matrix with the
#'   selection proportions of predictors.} \item{selectedY_full}{a binary matrix
#'   encoding stably selected outcomes.} \item{selpropY_full}{a matrix with the
#'   selection proportions of outcomes.} \item{coefX}{an array of loadings
#'   coefficients for the predictors (columns) over the fitted components (rows)
#'   and across the \code{K} resampling iterations (along the third dimension).}
#'   \item{coefY}{an array of loadings coefficients for the outcomes (columns)
#'   over the fitted components (rows) and across the \code{K} resampling
#'   iterations (along the third dimension). Only returned for supervised
#'   approaches (PLS).} \item{method}{a list with \code{type="bi_selection"},
#'   \code{implementation}, \code{family}, \code{resampling} and
#'   \code{PFER_method} values used for the run.} \item{params}{a list of input
#'   values used for the run. The datasets \code{xdata} and \code{ydata} are
#'   also included if \code{output_data=TRUE}.} The rows of \code{selectedX},
#'   \code{selectedY}, \code{selpropX} and \code{selpropY} correspond to the
#'   different combinations of parameters listed in \code{summary}. The rows of
#'   \code{selectedX_full}, \code{selectedY_full}, \code{selpropX_full} and
#'   \code{selpropY_full} correspond to the different combinations of parameters
#'   listed in \code{summary_full}.
#'
#' @family stability selection functions
#' @seealso \code{\link{SparsePLS}}, \code{\link{GroupPLS}},
#'   \code{\link{SparseGroupPLS}}, \code{\link{Resample}},
#'   \code{\link{StabilityScore}}
#'
#' @references \insertRef{sparsegroupPLS}{focus}
#'
#'   \insertRef{sparsePLS}{focus}
#'
#'   \insertRef{sparsePCA}{focus}
#'
#'   \insertRef{sparsePCASVD}{focus}
#'
#'   \insertRef{stabilityselectionMB}{focus}
#'
#'   \insertRef{stabilityselectionSS}{focus}
#'
#'   \insertRef{ourstabilityselection}{focus}
#'
#' @examples
#' \dontshow{
#'
#' # Data simulation
#' K <- 5
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#'
#' # sPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", K = K, ncomp = 2,
#'   LambdaX = 1:2,
#'   LambdaY = 1:2,
#'   implementation = SparsePLS
#' )
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateComponents(n = 50, pk = c(5, 5, 5))
#'
#' # sPCA: sparsity on X
#' stab <- BiSelection(
#'   xdata = simul$data,
#'   K = K, ncomp = 2,
#'   LambdaX = 1:2,
#'   implementation = SparsePCA
#' )
#' }
#'
#' \dontrun{
#' par(mar = c(12, 5, 1, 1))
#'
#'
#' ## Sparse Principal Component Analysis
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4))
#'
#' # sPCA: sparsity on X (unsupervised)
#' stab <- BiSelection(
#'   xdata = simul$data,
#'   ncomp = 3,
#'   LambdaX = 1:(ncol(simul$data) - 1),
#'   implementation = SparsePCA
#' )
#' print(stab)
#'
#' # Calibration plot
#' CalibrationPlot(stab)
#'
#' # Visualisation of the results
#' summary(stab)
#' plot(stab)
#' SelectedVariables(stab)
#'
#'
#' ## Sparse/Group Partial Least Squares
#'
#' # Data simulation (continuous outcomes)
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#' x <- simul$xdata
#' y <- simul$ydata
#'
#' # sPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#' plot(stab)
#'
#' # sPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   LambdaY = 1:(ncol(y) - 1),
#'   implementation = SparsePLS,
#'   n_cat = 2
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#' plot(stab)
#'
#' # sgPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y, K = 10,
#'   group_x = c(2, 8, 5),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.1),
#'   implementation = SparseGroupPLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#'
#' # sgPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y, K = 10,
#'   group_x = c(2, 8, 5), group_y = c(1, 2),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.2),
#'   LambdaY = 1:2, AlphaY = seq(0.1, 0.9, by = 0.2),
#'   implementation = SparseGroupPLS,
#'   n_cat = 2
#' )
#' CalibrationPlot(stab)
#' CalibrationPlot(stab,
#'   params = c("nx", "alphax", "ny", "alphay")
#' )
#' summary(stab)
#'
#' # gPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 8, 5),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2,
#'   implementation = GroupPLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#'
#' # gPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 8, 5), group_y = c(1, 2),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, LambdaY = 1:2,
#'   implementation = GroupPLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#'
#'
#' ## Sparse/Group PLS-DA (Discriminant Analysis)
#'
#' # Data simulation (categorical outcomes)
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = c(5, 5, 5), family = "binomial")
#' x <- simul$xdata
#' y <- simul$ydata
#'
#' # sPLS-DA: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y),
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#'
#' # sgPLS-DA: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y), K = 10,
#'   group_x = c(2, 8, 5),
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.1),
#'   implementation = SparseGroupPLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#'
#' # gPLS-DA: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y),
#'   group_x = c(2, 8, 5),
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:2,
#'   implementation = GroupPLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#' }
#' @export
BiSelection <- function(xdata, ydata = NULL, group_x = NULL, group_y = NULL,
                        LambdaX = NULL, LambdaY = NULL, AlphaX = NULL, AlphaY = NULL,
                        ncomp = 1, scale = TRUE,
                        pi_list = seq(0.6, 0.9, by = 0.01),
                        K = 100, tau = 0.5, seed = 1, n_cat = 3,
                        family = "gaussian", implementation = SparsePLS,
                        resampling = "subsampling", PFER_method = "MB",
                        PFER_thr = Inf, FDP_thr = Inf,
                        n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  # Defining Lambda if used with sparse PCA or PLS
  if (is.null(LambdaX)) {
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", "GroupPLS")) {
      LambdaX <- seq(1, length(group_x) - 1)
    }
    if (as.character(substitute(implementation)) %in% c("SparsePLS", "SparsePCA")) {
      LambdaX <- seq(1, ncol(xdata) - 1)
    }
  }

  # Preparing xdata
  xdata <- as.matrix(xdata)
  if (is.null(colnames(xdata))) {
    colnames(xdata) <- paste0("var", 1:ncol(xdata))
  }

  # Preparing ydata
  if (!is.null(ydata)) {
    if (is.vector(ydata) | is.factor(ydata)) {
      ydata <- matrix(ydata, ncol = 1)
    } else {
      if (family == "binomial") {
        ydata <- cbind(apply(ydata, 1, sum))
      }
    }
    if (is.null(colnames(ydata))) {
      colnames(ydata) <- paste0("outcome", 1:ncol(ydata))
    }
  }

  # Naming rows of xdata and ydata
  if (is.null(ydata)) {
    if (is.null(rownames(xdata))) {
      rownames(xdata) <- paste0("obs", 1:nrow(xdata))
    }
  } else {
    if (is.null(rownames(xdata)) & is.null(rownames(ydata))) {
      rownames(xdata) <- paste0("obs", 1:nrow(xdata))
      rownames(ydata) <- rownames(xdata)
    } else {
      if ((is.null(rownames(xdata))) & (!is.null(rownames(ydata)))) {
        rownames(xdata) <- rownames(ydata)
      }
      if ((!is.null(rownames(xdata))) & (is.null(rownames(ydata)))) {
        rownames(ydata) <- rownames(xdata)
      }
    }
  }

  if (as.character(substitute(implementation)) %in% c("SparsePLS", "GroupPLS")) {
    AlphaX <- AlphaY <- NULL
  }

  if (is.null(LambdaY)) {
    LambdaY <- NA
  }

  if (is.null(AlphaX)) {
    AlphaX <- NA
  }

  if (is.null(AlphaY)) {
    AlphaY <- NA
  }

  # Preparing empty objects to be filled
  selprop_x <- selprop_x_comp <- NULL
  selected_x <- selected_x_comp <- NULL
  selprop_y <- selprop_y_comp <- NULL
  selected_y <- selected_y_comp <- NULL
  params <- NULL
  params_comp <- matrix(NA, nrow = ncomp, ncol = 8)
  colnames(params_comp) <- c("comp", "nx", "alphax", "pix", "ny", "alphay", "piy", "S")

  for (comp in 1:ncomp) {
    # Preparing empty objects to be filled at current iteration (comp)
    tmp_selected_x <- matrix(NA, ncol = ncol(xdata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    tmp_selprop_x <- matrix(NA, ncol = ncol(xdata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    if (!is.null(ydata)) {
      if (family == "gaussian") {
        tmp_selected_y <- matrix(NA, ncol = ncol(ydata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
        tmp_selprop_y <- matrix(NA, ncol = ncol(ydata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
      } else {
        tmp_selected_y <- matrix(NA, ncol = length(unique(ydata)), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
        tmp_selprop_y <- matrix(NA, ncol = length(unique(ydata)), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
      }
    }
    tmp_params <- matrix(NA, ncol = 7, nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    colnames(tmp_params) <- c("nx", "alphax", "pix", "ny", "alphay", "piy", "S")

    # Initialisation of the run
    id <- 1
    if (verbose) {
      cat("\n")
      message(paste0("Component ", comp))
      if (as.character(substitute(implementation)) == "SparseGroupPLS") {
        pb <- utils::txtProgressBar(style = 3)
      }
    }

    # For loops over different grids of parameters (done internally in VariableSelection() for LambdaX)
    for (ny in LambdaY) {
      for (alphax in AlphaX) {
        for (alphay in AlphaY) {
          if (as.character(substitute(implementation)) == "SparsePCA") {
            if (family == "gaussian") {
              stab <- VariableSelection(
                xdata = xdata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                ncomp = comp, scale = scale, ...
              )
            } else {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                ncomp = comp, scale = scale, ...
              )
            }
          }
          if (as.character(substitute(implementation)) == "SparsePLS") {
            if (family == "gaussian") {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                keepY = NAToNULL(c(params_comp[1:comp, "ny"], ny)),
                ncomp = comp, scale = scale, ...
              )
            } else {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                ncomp = comp, scale = scale, ...
              )
            }
          }
          if (as.character(substitute(implementation)) == "SparseGroupPLS") {
            if (family == "gaussian") {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = FALSE,
                group_x = NAToNULL(group_x), group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                alpha.x = NAToNULL(c(params_comp[1:comp, "alphax"], alphax)),
                keepY = NAToNULL(c(params_comp[1:comp, "ny"], ny)),
                alpha.y = NAToNULL(c(params_comp[1:comp, "alphay"], alphay)),
                ncomp = comp, scale = scale, ...
              )
            } else {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = FALSE,
                group_x = NAToNULL(group_x), group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                alpha.x = NAToNULL(c(params_comp[1:comp, "alphax"], alphax)),
                ncomp = comp, scale = scale, ...
              )
            }
          }
          if (as.character(substitute(implementation)) == "GroupPLS") {
            if (family == "gaussian") {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                group_x = NAToNULL(group_x), group_y = group_y,
                group_penalisation = TRUE,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                keepY = NAToNULL(c(params_comp[1:comp, "ny"], ny)),
                ncomp = comp, scale = scale, ...
              )
            } else {
              stab <- VariableSelection(
                xdata = xdata, ydata = ydata,
                Lambda = LambdaX, pi_list = pi_list,
                K = K, tau = tau, seed = seed, n_cat = n_cat,
                family = family, implementation = implementation,
                resampling = resampling, PFER_method = PFER_method,
                PFER_thr = PFER_thr, FDP_thr = FDP_thr,
                n_cores = n_cores, output_data = FALSE, verbose = verbose,
                group_x = NAToNULL(group_x), group_y = group_y,
                group_penalisation = TRUE,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                ncomp = comp, scale = scale, ...
              )
            }
          }

          # Storing selections (X and Y)
          piy <- NULL
          if (!is.null(ydata)) {
            mycoefs <- Coefficients(stab, side = "Y", comp = comp)
          }
          for (i in 1:length(LambdaX)) {
            tmp_selprop_x[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- stab$selprop[i, ]
            tmp_selected_x[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- ifelse(stab$selprop[i, ] >= stab$P[i, ], yes = 1, no = 0)

            if (!is.null(ydata)) {
              tmpcoef <- mycoefs[i, , ]
              if (is.null(dim(tmpcoef))) {
                tmpcoef <- matrix(tmpcoef, nrow = dim(mycoefs)[2], ncol = dim(mycoefs)[3])
              }
              mytmp <- rep(NA, ifelse(length(dim(tmpcoef)) == 2, yes = nrow(tmpcoef), no = 1))
              for (l in 1:nrow(tmpcoef)) {
                mytmp[l] <- sum(tmpcoef[l, ] != 0) / length(tmpcoef[l, ])
              }
              tmp_selprop_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- mytmp
              if (any(mytmp != 1)) {
                # Computing stability score for Y variables
                if (as.character(substitute(implementation)) == "GroupPLS") {
                  hat_pi <- stab$params$pi_list[which.max(StabilityScore(mytmp,
                    pi_list = stab$params$pi_list, K = K,
                    group = NAToNULL(group_y)
                  ))]
                } else {
                  hat_pi <- stab$params$pi_list[which.max(StabilityScore(mytmp, pi_list = stab$params$pi_list, K = K))]
                }
                tmp_selected_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- ifelse(mytmp >= hat_pi, yes = 1, no = 0)
              } else {
                tmp_selected_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- 1
                hat_pi <- NA
              }
              piy <- c(piy, hat_pi)
            }
          }

          # Storing parameter values
          if (is.null(ydata)) {
            piy <- rep(NA, length(stab$Lambda))
          }
          tmp_params[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX)), ] <- cbind(stab$Lambda, alphax, stab$P, ny, alphay, piy, stab$S)

          # Incrementing loading bar and id
          if (verbose & (as.character(substitute(implementation)) == "SparseGroupPLS")) {
            utils::setTxtProgressBar(pb, id / (length(LambdaY) * length(AlphaX) * length(AlphaY)))
          }
          id <- id + 1
        }
      }
    }

    if (verbose & (as.character(substitute(implementation)) == "SparseGroupPLS")) {
      cat("\n")
    }

    # Filling big outputs
    params <- rbind(params, cbind(rep(comp, nrow(tmp_params)), tmp_params))
    selected_x <- rbind(selected_x, tmp_selected_x)
    selprop_x <- rbind(selprop_x, tmp_selprop_x)
    if (!is.null(ydata)) {
      selected_y <- rbind(selected_y, tmp_selected_y)
      selprop_y <- rbind(selprop_y, tmp_selprop_y)
    }

    # Filling best parameters by component
    params_comp[comp, "comp"] <- comp
    params_comp[comp, "nx"] <- tmp_params[which.max(tmp_params[, "S"]), "nx"]
    params_comp[comp, "alphax"] <- tmp_params[which.max(tmp_params[, "S"]), "alphax"]
    params_comp[comp, "pix"] <- tmp_params[which.max(tmp_params[, "S"]), "pix"]
    params_comp[comp, "ny"] <- tmp_params[which.max(tmp_params[, "S"]), "ny"]
    params_comp[comp, "alphay"] <- tmp_params[which.max(tmp_params[, "S"]), "alphay"]
    params_comp[comp, "piy"] <- tmp_params[which.max(tmp_params[, "S"]), "piy"]
    params_comp[comp, "S"] <- tmp_params[which.max(tmp_params[, "S"]), "S"]

    # Filling best selected/selection proportions by component
    selected_x_comp <- rbind(selected_x_comp, tmp_selected_x[which.max(tmp_params[, "S"]), ])
    selprop_x_comp <- rbind(selprop_x_comp, tmp_selprop_x[which.max(tmp_params[, "S"]), ])
    if (!is.null(ydata)) {
      selected_y_comp <- rbind(selected_y_comp, tmp_selected_y[which.max(tmp_params[, "S"]), ])
      selprop_y_comp <- rbind(selprop_y_comp, tmp_selprop_y[which.max(tmp_params[, "S"]), ])
    }
  }

  # Extracting X loadings coefficients
  coefs <- stab$Beta[ArgmaxId(stab)[1], , , drop = FALSE]
  coefs_X <- array(NA,
    dim = c(ncomp, sum(grepl("X_.*PC1", colnames(coefs))), dim(coefs)[3]),
    dimnames = list(
      paste0("PC", 1:ncomp),
      gsub("_PC.*", "", gsub(
        "X_", "",
        colnames(coefs)[grep("X_.*PC1", colnames(coefs))]
      )),
      paste0("iter", 1:dim(coefs)[3])
    )
  )
  for (i in 1:ncomp) {
    coefs_X[i, , ] <- coefs[1, grep(paste0("X_.*PC", i), colnames(coefs)), ]
  }

  # Extracting Y loadings coefficients
  if (!is.null(ydata)) {
    coefs_Y <- array(NA,
      dim = c(ncomp, sum(grepl("Y_.*PC1", colnames(coefs))), dim(coefs)[3]),
      dimnames = list(
        paste0("PC", 1:ncomp),
        gsub("_PC.*", "", gsub(
          "Y_", "",
          colnames(coefs)[grep("Y_.*PC1", colnames(coefs))]
        )),
        paste0("iter", 1:dim(coefs)[3])
      )
    )
    for (i in 1:ncomp) {
      coefs_Y[i, , ] <- coefs[1, grep(paste0("Y_.*PC", i), colnames(coefs)), ]
    }
  }

  # Excluding irrelevant columns
  colnames(params) <- c("comp", "nx", "alphax", "pix", "ny", "alphay", "piy", "S")
  params_comp <- params_comp[, which(!is.na(c(1, LambdaX[1], AlphaX[1], LambdaX[1], LambdaY[1], AlphaY[1], LambdaY[1], 1))), drop = FALSE]
  params <- params[, colnames(params_comp), drop = FALSE]

  # Assigning row and column names
  colnames(selected_x_comp) <- colnames(selprop_x_comp) <- colnames(selected_x) <- colnames(selprop_x) <- colnames(xdata)
  rownames(selected_x_comp) <- rownames(selprop_x_comp) <- paste0("comp", 1:ncomp)
  if (!is.null(ydata)) {
    if (family == "gaussian") {
      colnames(selected_y_comp) <- colnames(selprop_y_comp) <- colnames(selected_y) <- colnames(selprop_y) <- colnames(ydata)
    }
    rownames(selected_y_comp) <- rownames(selprop_y_comp) <- paste0("comp", 1:ncomp)
  }

  # Transposing selection status and proportion to be aligned with
  selected_x_comp <- t(selected_x_comp)
  selprop_x_comp <- t(selprop_x_comp)
  selected_x <- t(selected_x)
  selprop_x <- t(selprop_x)
  if (!is.null(ydata)) {
    selected_y_comp <- t(selected_y_comp)
    selprop_y_comp <- t(selprop_y_comp)
    if (family == "gaussian") {
      selected_y <- t(selected_y)
      selprop_y <- t(selprop_y)
    }
  }

  # Preparing outputs
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  if (!is.null(ydata)) {
    out <- list(
      summary = data.frame(params_comp),
      summary_full = data.frame(params),
      selectedX = selected_x_comp,
      selpropX = selprop_x_comp,
      selectedY = selected_y_comp,
      selpropY = selprop_y_comp,
      selectedX_full = selected_x,
      selpropX_full = selprop_x,
      selectedY_full = selected_y,
      selpropY_full = selprop_y,
      coefX = coefs_X,
      coefY = coefs_Y,
      methods = list(
        type = "bi_selection", implementation = as.character(substitute(implementation)),
        family = family, scale = scale,
        resampling = myresampling, PFER_method = PFER_method
      ),
      params = list(
        K = K, group_x = group_x, group_y = group_y,
        LambdaX = LambdaX, LambdaY = LambdaY,
        AlphaX = AlphaX, AlphaY = AlphaY,
        pi_list = pi_list,
        tau = tau, n_cat = n_cat, pk = ncol(xdata), n = nrow(xdata),
        PFER_thr = PFER_thr, FDP_thr = FDP_thr,
        seed = stab$seed
      )
    )
  } else {
    params_comp <- params_comp[, intersect(c("comp", "nx", "alphax", "pix", "S"), colnames(params_comp)), drop = FALSE]
    params <- params[, colnames(params_comp), drop = FALSE]
    out <- list(
      summary = data.frame(params_comp),
      summary_full = data.frame(params),
      selectedX = selected_x_comp,
      selpropX = selprop_x_comp,
      selectedX_full = selected_x,
      selpropX_full = selprop_x,
      coefX = coefs_X,
      methods = list(
        type = "bi_selection", implementation = as.character(substitute(implementation)),
        family = family, scale = scale,
        resampling = myresampling, PFER_method = PFER_method
      ),
      params = list(
        K = K,
        LambdaX = LambdaX,
        AlphaX = AlphaX,
        pi_list = pi_list,
        tau = tau, n_cat = n_cat, pk = ncol(xdata), n = nrow(xdata),
        PFER_thr = PFER_thr, FDP_thr = FDP_thr,
        seed = stab$seed
      )
    )
  }

  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata, ydata = ydata))
  }

  # Defining the class
  class(out) <- "bi_selection"

  return(out)
}

#' Variable selection for predictors and outcomes
#'
#' Runs stability selection regression models with different combinations of
#' parameters controlling the sparsity in PLS models and thresholds in selection
#' proportions. These parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features).
#'
#' @inheritParams VariableSelection
#' @param group_x vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group. Only used with
#'   \code{implementation="SparseGroupPLS"}.
#' @param group_y optional vector encoding the grouping structure among
#'   outcomes. This argument indicates the number of variables in each group.
#'   Only used with \code{implementation="SparseGroupPLS"}.
#' @param LambdaX matrix of parameters controlling the number of selected
#'   variables (sparse PLS) or groups (sparse group PLS) in X.
#' @param LambdaY matrix of parameters controlling the number of selected
#'   variables (sparse PLS) or groups (sparse group PLS) in Y. Only used with
#'   \code{family="gaussian"}.
#' @param AlphaX matrix of parameters controlling the level of sparsity within
#'   groups (sparse group PLS) in X. Only used with
#'   \code{implementation="SparseGroupPLS"}.
#' @param AlphaY matrix of parameters controlling the level of sparsity within
#'   groups (sparse group PLS) in X. Only used with
#'   \code{implementation="SparseGroupPLS"} and \code{family="gaussian"}.
#' @param ncomp number of components.
#'
#' @family stability selection functions
#' @seealso \code{\link{SparsePLS}}, \code{\link{GroupPLS}}, \code{\link{SparseGroupPLS}}
#'
#' @examples
#' \dontshow{
#'
#' # Data simulation
#' K <- 5
#' pk <- 15
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
#' ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
#' colnames(ydata) <- paste0("outcome", 1:4)
#' x <- simul$X
#' y <- ydata
#'
#' # sPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", K = K, ncomp = 2,
#'   LambdaX = 1:2,
#'   LambdaY = 1:2,
#'   implementation = SparsePLS
#' )
#' }
#'
#' \dontrun{
#'
#' # Data simulation (continuous outcomes)
#' pk <- 15
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
#' ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
#' colnames(ydata) <- paste0("outcome", 1:4)
#' x <- simul$X
#' y <- ydata
#'
#' # sPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#'
#' # sPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   LambdaY = 1:(ncol(y) - 1),
#'   implementation = SparsePLS
#' )
#'
#' # sgPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 10, 8),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.1),
#'   implementation = SparseGroupPLS
#' )
#'
#' # sgPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 10, 8), group_y = c(1, 3),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.1),
#'   LambdaY = 1:2, AlphaY = seq(0.1, 0.9, by = 0.1),
#'   implementation = SparseGroupPLS
#' )
#'
#' # gPLS: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 10, 8),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2,
#'   implementation = GroupPLS
#' )
#'
#' # gPLS: sparsity on both X and Y
#' stab <- BiSelection(
#'   xdata = x, ydata = y,
#'   group_x = c(2, 10, 8), group_y = c(1, 3),
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:2, LambdaY = 1:2,
#'   implementation = GroupPLS
#' )
#'
#' # Data simulation (categorical outcomes)
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#' x <- simul$X
#' y <- cbind(simul$Y, matrix(sample(c(0, 1), size = 200 * 3, replace = TRUE), ncol = 3))
#' y <- apply(y, 1, sum)
#'
#' # sPLS-DA: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y),
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   implementation = SparsePLS
#' )
#'
#' # sgPLS-DA: sparsity on X
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y),
#'   group_x = c(2, 10, 8), K = 10,
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:2, AlphaX = seq(0.1, 0.9, by = 0.1),
#'   implementation = SparseGroupPLS
#' )
#'
#' # gPLS-DA: sparsity on Y
#' stab <- BiSelection(
#'   xdata = x, ydata = cbind(y),
#'   group_x = c(2, 10, 8),
#'   family = "binomial", ncomp = 3,
#'   LambdaX = 1:2,
#'   implementation = GroupPLS
#' )
#' }
#' @export
BiSelection <- function(xdata, ydata, group_x = NULL, group_y = NULL,
                        LambdaX = NULL, LambdaY = NULL, AlphaX = NULL, AlphaY = NULL, ncomp = 1,
                        pi_list = seq(0.6, 0.9, by = 0.01),
                        K = 100, tau = 0.5, seed = 1, n_cat = 3,
                        family = "gaussian", implementation = SparsePLS,
                        resampling = "subsampling", PFER_method = "MB",
                        PFER_thr = Inf, FDP_thr = Inf,
                        n_cores = 1, output_data = FALSE, verbose = TRUE, ...) {
  if (is.null(LambdaX)) {
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", "GroupPLS")) {
      LambdaX <- 1:length(group_x)
    }
    if (as.character(substitute(implementation)) == "SparsePLS") {
      LambdaX <- 1:ncol(xdata)
    }
  }

  # Reformatting inputs
  if (family == "binomial") {
    if (is.vector(ydata)) {
      ydata <- cbind(ydata)
    } else {
      ydata <- apply(ydata, 1, sum)
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
  colnames(params_comp) <- c("comp", "nx", "alphax", "pix", "ny", "alphay", "piy", "stability_score")

  for (comp in 1:ncomp) {
    # Preparing empty objects to be filled at current iteration (comp)
    tmp_selected_x <- matrix(NA, ncol = ncol(xdata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    tmp_selprop_x <- matrix(NA, ncol = ncol(xdata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    if (family == "gaussian") {
      tmp_selected_y <- matrix(NA, ncol = ncol(ydata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
      tmp_selprop_y <- matrix(NA, ncol = ncol(ydata), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    } else {
      tmp_selected_y <- matrix(NA, ncol = length(unique(ydata)), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
      tmp_selprop_y <- matrix(NA, ncol = length(unique(ydata)), nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    }
    tmp_params <- matrix(NA, ncol = 7, nrow = length(LambdaX) * length(LambdaY) * length(AlphaX) * length(AlphaY))
    colnames(tmp_params) <- c("nx", "alphax", "pix", "ny", "alphay", "piy", "stability_score")

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
                ncomp = comp, ...
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
                ncomp = comp, ...
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
                group_x = group_x, group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                alpha.x = NAToNULL(c(params_comp[1:comp, "alphax"], alphax)),
                keepY = NAToNULL(c(params_comp[1:comp, "ny"], ny)),
                alpha.y = NAToNULL(c(params_comp[1:comp, "alphay"], alphay)),
                ncomp = comp, ...
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
                group_x = group_x, group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                alpha.x = NAToNULL(c(params_comp[1:comp, "alphax"], alphax)),
                ncomp = comp, ...
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
                group_x = group_x, group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                keepY = NAToNULL(c(params_comp[1:comp, "ny"], ny)),
                ncomp = comp, ...
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
                group_x = group_x, group_y = group_y,
                keepX_previous = NAToNULL(params_comp[1:comp, "nx"]),
                ncomp = comp, ...
              )
            }
          }

          # Storing selections (X and Y)
          piy <- NULL
          for (i in 1:length(LambdaX)) {
            tmp_selprop_x[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- stab$selprop[i, ]
            tmp_selected_x[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- ifelse(stab$selprop[i, ] >= stab$P[i, ], yes = 1, no = 0)
            mytmp <- apply(Coefficients(stab, side = "Y", comp = comp)[i, , ], 1, FUN = function(z) {
              sum(z != 0) / length(z)
            })
            tmp_selprop_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- mytmp
            if (any(mytmp != 1)) {
              hat_pi <- stab$params$pi_list[which.max(StabilityScore(mytmp, pi_list = stab$params$pi_list, K = K))]
              tmp_selected_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- ifelse(mytmp >= hat_pi, yes = 1, no = 0)
            } else {
              tmp_selected_y[seq((id - 1) * length(LambdaX) + 1, id * length(LambdaX))[i], ] <- 1
              hat_pi <- NA
            }
            piy <- c(piy, hat_pi)
          }

          # Storing parameter values
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
    selected_y <- rbind(selected_y, tmp_selected_y)
    selprop_y <- rbind(selprop_y, tmp_selprop_y)

    # Filling best parameters by component
    params_comp[comp, "comp"] <- comp
    params_comp[comp, "nx"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "nx"]
    params_comp[comp, "alphax"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "alphax"]
    params_comp[comp, "pix"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "pix"]
    params_comp[comp, "ny"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "ny"]
    params_comp[comp, "alphay"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "alphay"]
    params_comp[comp, "piy"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "piy"]
    params_comp[comp, "stability_score"] <- tmp_params[which.max(tmp_params[, "stability_score"]), "stability_score"]

    # Filling best selected/selection proportions by component
    selected_x_comp <- rbind(selected_x_comp, tmp_selected_x[which.max(tmp_params[, "stability_score"]), ])
    selprop_x_comp <- rbind(selprop_x_comp, tmp_selprop_x[which.max(tmp_params[, "stability_score"]), ])
    selected_y_comp <- rbind(selected_y_comp, tmp_selected_y[which.max(tmp_params[, "stability_score"]), ])
    selprop_y_comp <- rbind(selprop_y_comp, tmp_selprop_y[which.max(tmp_params[, "stability_score"]), ])
  }
  colnames(params) <- c("comp", "nx", "alphax", "pix", "ny", "alphay", "piy", "stability_score")

  # Assigning column names
  colnames(selected_x_comp) <- colnames(selprop_x_comp) <- colnames(selected_x) <- colnames(selprop_x) <- colnames(xdata)
  colnames(selected_y_comp) <- colnames(selprop_y_comp) <- colnames(selected_y) <- colnames(selprop_y) <- colnames(ydata)

  # Preparing outputs
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  out <- list(
    summary = params_comp,
    summary_full = params,
    selectedX = selected_x_comp,
    selpropX = selprop_x_comp,
    selectedY = selected_y_comp,
    selpropY = selprop_y_comp,
    selectedX_full = selected_x,
    selpropX_full = selprop_x,
    selectedY_full = selected_y,
    selpropY_full = selprop_y,
    methods = list(
      implementation = as.character(substitute(implementation)), family = family,
      resampling = myresampling, PFER_method = PFER_method
    ),
    params = list(
      K = K, group_x = group_x, group_y = group_y,
      LambdaX = LambdaX, LambdaY = LambdaY,
      AlphaX = AlphaX, AlphaY = AlphaY,
      pi_list = pi_list,
      tau = tau, n_cat = n_cat, pk = ncol(xdata), n = nrow(xdata),
      PFER_thr = PFER_thr, FDP_thr = FDP_thr,
      seed = seed
    )
  )

  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata, ydata = ydata))
  }

  return(out)
}

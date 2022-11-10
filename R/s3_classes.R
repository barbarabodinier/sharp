#' @export
print.variable_selection <- function(x, ...) {
  cat(paste0("Stability selection using function ", x$methods$implementation, " with ", x$methods$family, " family."))
  cat("\n")
  if (x$methods$resampling == "subsampling") {
    cat(paste0("The model was run using ", x$params$K, " subsamples of ", x$params$tau * 100, "% of the observations."))
  } else {
    cat(paste0("The model was run using ", x$params$K, " bootstrap samples."))
  }
  cat("\n")
}


#' @export
print.graphical_model <- function(x, ...) {
  cat(paste0("Stability selection using function ", x$methods$implementation, "."))
  cat("\n")
  if (x$methods$resampling == "subsampling") {
    cat(paste0("The model was run using ", x$params$K, " subsamples of ", x$params$tau * 100, "% of the observations."))
  } else {
    cat(paste0("The model was run using ", x$params$K, " bootstrap samples."))
  }
  cat("\n")
}


#' @export
print.bi_selection <- function(x, ...) {
  cat(paste0("Stability selection using function ", x$methods$implementation, " with ", x$methods$family, " family."))
  cat("\n")
  if (x$methods$resampling == "subsampling") {
    cat(paste0("The model was run using ", x$params$K, " subsamples of ", x$params$tau * 100, "% of the observations."))
  } else {
    cat(paste0("The model was run using ", x$params$K, " bootstrap samples."))
  }
  cat("\n")
}


#' @export
summary.variable_selection <- function(object, ...) {
  cat(paste0(
    "Calibrated parameters: lambda = ",
    formatC(Argmax(object)[1, 1], format = "f", digits = 3),
    " and pi = ",
    formatC(Argmax(object)[1, 2], format = "f", digits = 3)
  ))
  cat("\n")
  cat("\n")
  cat(paste0(
    "Maximum stability score: ",
    formatC(max(object$S, na.rm = TRUE), format = "f", digits = 3)
  ))
  cat("\n")
  cat("\n")
  cat(paste0(
    "Number of selected variable(s): ",
    sum(SelectedVariables(object))
  ))
  cat("\n")
}


#' @export
summary.graphical_model <- function(object, ...) {
  if (ncol(object$S) > 1) {
    cat(paste0("Calibrated parameters:"))
    cat("\n")
    for (k in 1:ncol(object$S)) {
      cat(paste0(
        "Block ", k, ": lambda = ",
        formatC(Argmax(object)[k, 1], format = "f", digits = 3),
        " and pi = ",
        formatC(Argmax(object)[k, 2], format = "f", digits = 3)
      ))
      cat("\n")
    }
    cat("\n")
    cat("Maximum stability scores: ")
    cat("\n")
    for (k in 1:ncol(object$S)) {
      cat(paste0(
        "Block ", k, ": ",
        formatC(max(object$S[, k], na.rm = TRUE), format = "f", digits = 3)
      ))
      cat("\n")
    }
    cat("\n")
    cat("Number of selected edge(s): ")
    cat("\n")
    adjacency <- Adjacency(object)
    adjacency <- adjacency[upper.tri(adjacency)]
    bigblocks <- BlockMatrix(pk = object$params$pk)
    bigblocks <- bigblocks[upper.tri(bigblocks)]
    for (k in 1:ncol(object$S)) {
      cat(paste0(
        "Block ", k, ": ",
        round(sum(adjacency[bigblocks == k]))
      ))
      cat("\n")
    }
    cat(paste0(
      "Total: ",
      sum(Adjacency(object)) / 2
    ))
  } else {
    cat(paste0(
      "Calibrated parameters: lambda = ",
      formatC(Argmax(object)[1, 1], format = "f", digits = 3),
      " and pi = ",
      formatC(Argmax(object)[1, 2], format = "f", digits = 3)
    ))
    cat("\n")
    cat("\n")
    cat(paste0(
      "Maximum stability score: ",
      formatC(max(object$S[, 1], na.rm = TRUE), format = "f", digits = 3)
    ))
    cat("\n")
    cat("\n")
    cat(paste0(
      "Number of selected edge(s): ",
      sum(Adjacency(object)) / 2
    ))
  }
  cat("\n")
}


#' @export
summary.bi_selection <- function(object, ...) {
  cat(paste0("Calibrated parameters (X):"))
  cat("\n")
  for (k in 1:nrow(object$summary)) {
    if ("alphax" %in% colnames(object$summary)) {
      cat(paste0(
        "Component ", k, ": n = ",
        formatC(object$summary[k, "nx"], format = "f", digits = 3),
        ", alpha = ",
        formatC(object$summary[k, "alphax"], format = "f", digits = 3),
        " and pi = ",
        formatC(object$summary[k, "pix"], format = "f", digits = 3)
      ))
    } else {
      cat(paste0(
        "Component ", k, ": n = ",
        formatC(object$summary[k, "nx"], format = "f", digits = 3),
        " and pi = ",
        formatC(object$summary[k, "pix"], format = "f", digits = 3)
      ))
    }
    cat("\n")
  }
  if ("ny" %in% colnames(object$summary)) {
    cat("\n")
    cat(paste0("Calibrated parameters (Y):"))
    cat("\n")
    for (k in 1:nrow(object$summary)) {
      if ("alphay" %in% colnames(object$summary)) {
        cat(paste0(
          "Component ", k, ": n = ",
          formatC(object$summary[k, "ny"], format = "f", digits = 3),
          ", alpha = ",
          formatC(object$summary[k, "alphay"], format = "f", digits = 3),
          " and pi = ",
          formatC(object$summary[k, "piy"], format = "f", digits = 3)
        ))
      } else {
        cat(paste0(
          "Component ", k, ": n = ",
          formatC(object$summary[k, "ny"], format = "f", digits = 3),
          " and pi = ",
          formatC(object$summary[k, "piy"], format = "f", digits = 3)
        ))
      }
      cat("\n")
    }
  }
  cat("\n")
  if (nrow(object$summary) > 1) {
    cat("Maximum stability scores (X): ")
  } else {
    cat("Maximum stability score (X): ")
  }
  cat("\n")
  for (k in 1:nrow(object$summary)) {
    cat(paste0(
      "Component ", k, ": ",
      formatC(max(object$summary[k, "S"], na.rm = TRUE), format = "f", digits = 3)
    ))
    cat("\n")
  }
  cat("\n")
  cat("Number of selected variable(s) (X): ")
  cat("\n")
  for (k in 1:nrow(object$summary)) {
    cat(paste0(
      "Component ", k, ": ",
      round(sum(object$selectedX[, k]))
    ))
    cat("\n")
  }
  if ("ny" %in% colnames(object$summary)) {
    cat("\n")
    cat("Number of selected variable(s) (Y): ")
    cat("\n")
    for (k in 1:nrow(object$summary)) {
      cat(paste0(
        "Component ", k, ": ",
        round(sum(object$selectedY[, k]))
      ))
      cat("\n")
    }
  }
  cat("\n")
}


#' @export
plot.variable_selection <- function(x, ...) {
  selprop <- SelectionProportions(x)
  selprop <- sort(selprop, decreasing = TRUE)
  plot(selprop,
    type = "h",
    xlab = "", ylab = "Selection proportion",
    las = 1, ylim = c(0, 1), xaxt = "n",
    col = ifelse(SelectedVariables(x)[names(selprop)],
      yes = "navy", no = "grey"
    )
  )
  graphics::axis(
    side = 1, at = 1:length(selprop),
    labels = names(selprop), las = 2
  )
  graphics::abline(h = Argmax(x)[1, 2], lty = 2, col = "darkred")
}


#' @export
plot.graphical_model <- function(x, ...) {
  igraph::plot.igraph(Graph(x), ...)
}


#' @export
plot.bi_selection <- function(x, ...) {
  igraph::plot.igraph(Graph(x), ...)
}


#' @export
coef.variable_selection <- function(object, ...) {
  # Checking inputs
  if (!object$methods$family %in% c("gaussian", "binomial", "multinomial", "mgaussian", "cox")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian', 'binomial', 'multinomial', 'mgaussian' or 'cox'.")
  }

  # Extracting index of calibrated parameter
  argmax_id <- ArgmaxId(stability = object)[1]

  # Extracting beta coefficients
  if (object$methods$family %in% c("gaussian", "binomial", "cox")) {
    beta <- t(object$Beta[argmax_id, , ])
    rownames(beta) <- paste0("iter", 1:nrow(beta))
  }
  if (object$methods$family %in% c("multinomial", "mgaussian")) {
    tmpbeta <- object$Beta[argmax_id, , , ]
    beta <- array(NA, dim = c(dim(tmpbeta)[2], dim(tmpbeta)[1], dim(tmpbeta)[3]))
    for (k in 1:dim(tmpbeta)[3]) {
      beta[, , k] <- t(tmpbeta[, , k])
    }
    dimnames(beta) <- list(
      paste0("iter", 1:nrow(beta)),
      dimnames(object$Beta)[[2]],
      dimnames(object$Beta)[[4]]
    )
  }
  # Intercept is not included but could be obtained from Ensemble() for "gaussian" or "binomial"

  return(beta)
}


#' @export
predict.variable_selection <- function(object, xdata, ydata, newdata = NULL, method = c("ensemble", "refit"), ...) {
  # Checking inputs
  if (!object$methods$family %in% c("gaussian", "binomial", "multinomial", "cox")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian', 'binomial', 'multinomial' or 'cox'.")
  } else {
    if (method[[1]] == "ensemble") {
      if (!object$methods$family %in% c("gaussian", "binomial")) {
        method <- "refit"
        message("Predictions from ensemble models is only available for the following families for regression models: 'gaussian' or 'binomial'. Predicted values are obtained from refitting.")
      }
    }
  }

  # Using the same data if not provided
  if (is.null(newdata)) {
    newdata <- xdata
  }

  # Predictions from ensemble model
  if (method[1] == "ensemble") {
    ensemble <- Ensemble(
      stability = object,
      xdata = xdata,
      ydata = ydata
    )
    yhat <- EnsemblePredictions(
      ensemble = ensemble,
      xdata = newdata,
      ...
    )
  }

  # Predictions from refitted model
  if (method[1] == "refit") {
    refitted <- Refit(xdata = xdata, ydata = ydata, stability = object)
    yhat <- stats::predict(object = refitted, newdata = as.data.frame(xdata), ...)
    yhat <- cbind(yhat)
  }
  return(yhat)
}


#' Receiver Operating Characteristic (ROC) band
#'
#' Plots the True Positive Rate (TPR) as a function of the False Positive Rate
#' (FPR) for different thresholds in predicted probabilities. If the results
#' from multiple ROC analyses are provided (e.g. output of
#' \code{\link{ExplanatoryPerformance}} with large \code{K}), the point-wise
#' median is represented and flanked by a transparent band defined by point-wise
#' \code{quantiles}.
#'
#' @param x output of \code{\link{ExplanatoryPerformance}}.
#' @param col_band colour of the band defined by point-wise \code{quantiles}.
#' @param alpha level of opacity for the band.
#' @param quantiles point-wise quantiles of the performances defining the band.
#' @param add logical indicating if the curve should be added to the current
#'   plot.
#' @param ... additional plotting arguments (see \code{\link[graphics]{par}}).
#'
#' @return A base plot.
#'
#' @seealso \code{\link{ExplanatoryPerformance}}
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 1000, pk = 10,
#'   family = "binomial", ev_xy = 0.7
#' )
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Evaluation of the performances on refitted models (K=1)
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, n_thr = NULL
#' )
#' plot(roc)
#'
#' # Using more refitting/test splits
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 100
#' )
#' boxplot(roc$AUC, ylab = "AUC")
#' plot(roc)
#' }
#' @export
plot.roc_band <- function(x,
                          col_band = NULL,
                          alpha = 0.5,
                          quantiles = c(0.05, 0.95),
                          add = FALSE,
                          ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Defining default parameters if not provided
  if (!"xlim" %in% names(extra_args)) {
    extra_args$xlim <- c(0, 1)
  }
  if (!"ylim" %in% names(extra_args)) {
    extra_args$ylim <- c(0, 1)
  }
  if (!"lwd" %in% names(extra_args)) {
    extra_args$lwd <- 2
  }
  if (!"xlab" %in% names(extra_args)) {
    extra_args$xlab <- "False Positive Rate"
  }
  if (!"ylab" %in% names(extra_args)) {
    extra_args$ylab <- "True Positive Rate"
  }
  if (!"las" %in% names(extra_args)) {
    extra_args$las <- 1
  }
  if (!"cex.lab" %in% names(extra_args)) {
    extra_args$cex.lab <- 1.3
  }
  if (!"col" %in% names(extra_args)) {
    extra_args$col <- "red"
  }

  # Extracting the number of iterations
  niter <- length(x$AUC)

  # Defining the band colour
  if (is.null(col_band)) {
    col_band <- extra_args$col
  }

  # Initialising the plot
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = plot)
  if (!add) {
    do.call(plot, args = c(list(x = NULL), tmp_extra_args))
    graphics::abline(0, 1, lty = 3)
  }

  # Defining quantile bands
  if (nrow(x$FPR) > 1) {
    xseq <- apply(x$FPR, 2, FUN = function(x) {
      sort(x)[rev(quantiles) * niter]
    })
    yseq <- apply(x$TPR, 2, FUN = function(x) {
      sort(x)[quantiles * niter]
    })
    graphics::polygon(c(xseq[1, ], rev(xseq[2, ])),
      c(yseq[1, ], rev(yseq[2, ])),
      col = grDevices::adjustcolor(col = col_band, alpha.f = alpha),
      border = NA
    )
  }

  # Adding the point-wise average
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = graphics::lines)
  do.call(graphics::lines, args = c(
    list(
      x = apply(x$FPR, 2, stats::median),
      y = apply(x$TPR, 2, stats::median)
    ),
    tmp_extra_args
  ))
}


#' @export
summary.incremental <- function(object, ...) {
  cat(paste0("Performances of recalibrated models:"))
  cat("\n")
  cat("\n")
  mat <- PlotIncremental(object, output_data = TRUE, ...)
  for (i in 1:ncol(mat)) {
    cat(paste0(
      ifelse(i == 1, yes = "  ", no = "+ "),
      colnames(mat)[i],
      ": ",
      formatC(mat[1, i], format = "f", digits = 3)
    ))
    cat("\n")
  }
}


#' @export
plot.incremental <- function(x, ...) {
  PlotIncremental(x, ...)
}

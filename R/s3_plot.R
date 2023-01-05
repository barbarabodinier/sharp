#' Plot of selection proportions
#'
#' Makes a barplot of selection proportions in decreasing order. See examples in
#' \code{\link{VariableSelection}}.
#'
#' @param x output of \code{\link{VariableSelection}}.
#' @param col vector of colours by stable selection status.
#' @param col.axis optional vector of label colours by stable selection status.
#' @param col.thr threshold colour.
#' @param lty.thr threshold line type as \code{lty} in
#'   \code{\link[graphics]{par}}.
#' @param n_predictors number of predictors to display.
#' @param ... additional plotting arguments (see \code{\link[graphics]{par}}).
#'
#' @return A plot.
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @export
plot.variable_selection <- function(x,
                                    col = c("red", "grey"),
                                    col.axis = NULL,
                                    col.thr = "darkred",
                                    lty.thr = 2,
                                    n_predictors = NULL,
                                    ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Checking inputs
  if (length(col) != 2) {
    col <- rep(col[1], 2)
  }
  if (is.null(col.axis)) {
    col.axis <- col
  }

  # Defining default settings
  if (!"type" %in% names(extra_args)) {
    extra_args$type <- "h"
  }
  if (!"xlab" %in% names(extra_args)) {
    extra_args$xlab <- ""
  }
  if (!"ylab" %in% names(extra_args)) {
    extra_args$ylab <- "Selection proportion"
  }
  if (!"las" %in% names(extra_args)) {
    extra_args$las <- 2
  }
  if (!"ylim" %in% names(extra_args)) {
    extra_args$ylim <- c(0, 1)
  }
  if (!"cex.lab" %in% names(extra_args)) {
    extra_args$cex.lab <- 1.5
  }
  if (!"cex.axis" %in% names(extra_args)) {
    extra_args$cex.axis <- 1
  }

  # Extracting selection proportions
  selprop <- SelectionProportions(x)

  # Defining colours
  mycolours <- ifelse(SelectedVariables(x) == 1, yes = col[1], no = col[2])
  mycolours_axis <- mycolours

  # Re-ordering by decreasing selection proportion
  ids <- sort.list(selprop, decreasing = TRUE)
  if (is.null(n_predictors)) {
    n_predictors <- length(ids)
  }
  n_predictors <- min(n_predictors, length(ids))
  ids <- ids[1:n_predictors]
  selprop <- selprop[ids]
  mycolours <- mycolours[ids]
  mycolours_axis <- mycolours_axis[ids]

  # Extracting relevant extra arguments
  tmp_extra_args <- extra_args
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("xaxt", "col")]

  # Making plot
  do.call(base::plot, args = c(
    list(
      x = selprop,
      xaxt = "n",
      col = mycolours
    ),
    tmp_extra_args
  ))

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = graphics::axis)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("side", "at", "labels", "col.axis", "col.ticks", "type")]

  # Adding axis
  for (i in 1:length(selprop)) {
    do.call(graphics::axis, args = c(
      list(
        side = 1,
        at = i,
        labels = names(selprop)[i],
        col.axis = mycolours_axis[i],
        col.ticks = mycolours_axis[i]
      ),
      tmp_extra_args
    ))
  }
  graphics::abline(h = Argmax(x)[1, 2], lty = lty.thr, col = col.thr)
}


#' @export
plot.graphical_model <- function(x, ...) {
  igraph::plot.igraph(Graph(x), ...)
}


#' @export
plot.bi_selection <- function(x, ...) {
  igraph::plot.igraph(Graph(x), ...)
}


#' Consensus matrix heatmap
#'
#' Creates a heatmap of the (calibrated) consensus matrix. See examples in
#' \code{\link{Clustering}}.
#'
#' @inheritParams Clusters
#' @param x output of \code{\link{Clustering}}.
#' @param theta optional vector of cluster membership. If provided, the ordering
#'   of the items should be the same as in \code{\link{Clusters}}. This argument
#'   is used to re-order the consensus matrix.
#' @param theta_star optional vector of true cluster membership. If provided,
#'   the ordering of the items should be the same as in \code{\link{Clusters}}.
#'   This argument is used to define item colours.
#' @param lines logical indicating if lines separating the clusters provided in
#'   \code{theta} should be displayed.
#' @param col.lines colour of the lines separating the clusters.
#' @param lwd.lines width of the lines separating the clusters.
#' @param tick logical indicating if axis tickmarks should be displayed.
#' @param axes logical indicating if item labels should be displayed.
#' @param col.axis optional vector of cluster colours.
#' @param cex.axis font size for axes.
#' @param xlas orientation of labels on the x-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylas orientation of labels on the y-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param bty character string indicating if the box around the plot should be
#'   drawn. Possible values include: \code{"o"} (default, the box is drawn), or
#'   \code{"n"} (no box).
#' @param ... additional arguments passed to \code{\link{Heatmap}}.
#'
#' @return A heatmap.
#'
#' @export
plot.clustering <- function(x,
                            linkage = "complete",
                            argmax_id = NULL,
                            theta = NULL,
                            theta_star = NULL,
                            lines = TRUE,
                            col.lines = c("blue"),
                            lwd.lines = 2,
                            tick = TRUE,
                            axes = TRUE,
                            col.axis = NULL,
                            cex.axis = 1,
                            xlas = 2,
                            ylas = 2,
                            bty = "n",
                            ...) {
  # Defining theta if not provided
  if (is.null(theta)) {
    if (is.null(argmax_id)) {
      argmax_id <- ArgmaxId(stability = x)
    }
    coprop <- ConsensusMatrix(stability = x, argmax_id = argmax_id[1])
    shclust <- stats::hclust(stats::as.dist(1 - coprop), method = linkage)
    theta <- stats::cutree(shclust, k = ceiling(x$nc[argmax_id[1], 1]))
    ids <- shclust$order
  } else {
    theta <- as.numeric(as.factor(theta))
    ids <- sort.list(theta)
  }

  # Defining theta_star if not provided
  if (is.null(theta_star)) {
    theta_star <- theta
  } else {
    theta_star <- as.numeric(as.factor(theta_star))
  }

  # Ordering by theta
  mat <- ConsensusMatrix(stability = x, argmax_id = argmax_id)
  theta <- theta[ids]
  theta_star <- theta_star[ids]
  mat <- mat[ids, ids]

  # Preparing heatmap
  Heatmap(mat = mat, axes = FALSE, bty = bty, ...)

  # Adding separation lines based on theta
  if (lines) {
    thr <- which(!duplicated(theta))
    graphics::abline(h = length(theta) - thr[-1] + 1, col = col.lines, lwd = lwd.lines)
    graphics::abline(v = thr[-1] - 1, col = col.lines, lwd = lwd.lines)
  }

  # Defining axis colours if not provided
  if (is.null(col.axis)) {
    col.axis <- randomcoloR::distinctColorPalette(k = length(unique(theta_star)))
  }

  # Adding axes and ticks
  for (k in 1:nrow(mat)) {
    graphics::axis(
      side = 2,
      at = k - 0.5,
      labels = ifelse(axes, yes = rev(rownames(mat))[k], no = ""),
      las = ylas,
      col.axis = col.axis[rev(theta_star)[k]],
      col.ticks = col.axis[rev(theta_star)[k]],
      cex.axis = cex.axis, tick = tick
    )
    graphics::axis(
      side = 1,
      at = k - 0.5,
      labels = ifelse(axes, yes = colnames(mat)[k], no = ""),
      las = xlas,
      col.axis = col.axis[theta_star[k]],
      col.ticks = col.axis[theta_star[k]],
      cex.axis = cex.axis, tick = tick
    )
  }
}


#' Receiver Operating Characteristic (ROC) band
#'
#' Plots the True Positive Rate (TPR) as a function of the False Positive Rate
#' (FPR) for different thresholds in predicted probabilities. If the results
#' from multiple ROC analyses are provided (e.g. output of
#' \code{\link{ExplanatoryPerformance}} with large \code{K}), the point-wise
#' median is represented and flanked by a transparent band defined by point-wise
#' \code{quantiles}. See examples in \code{\link{ROC}} and
#' \code{\link{ExplanatoryPerformance}}.
#'
#' @param x output of \code{\link{ROC}} or \code{\link{ExplanatoryPerformance}}.
#' @param col_band colour of the band defined by point-wise \code{quantiles}.
#' @param alpha level of opacity for the band.
#' @param quantiles point-wise quantiles of the performances defining the band.
#' @param add logical indicating if the curve should be added to the current
#'   plot.
#' @param ... additional plotting arguments (see \code{\link[graphics]{par}}).
#'
#' @return A base plot.
#'
#' @seealso \code{\link{ROC}}, \code{\link{ExplanatoryPerformance}}
#'
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


#' Visualisation of incremental performance
#'
#' Represents prediction performances upon sequential inclusion of the
#' predictors in a logistic or Cox regression model as produced by
#' \code{\link{Incremental}}. The median and \code{quantiles} of the performance
#' metric are reported. See examples in \code{\link{Incremental}}.
#'
#' @inheritParams CalibrationPlot
#' @param x output of \code{\link{Incremental}}.
#' @param quantiles quantiles defining the lower and upper bounds.
#' @param sfrac size of the end bars, as in \code{\link[plotrix]{plotCI}}.
#' @param col vector of colours by stable selection status.
#' @param col.axis optional vector of label colours by stable selection status.
#' @param xcex.axis size of labels along the x-axis.
#' @param ycex.axis size of labels along the y-axis.
#' @param output_data logical indicating if the median and quantiles should be
#'   returned in a matrix.
#' @param ... additional plotting arguments (see \code{\link[graphics]{par}}).
#'
#' @return A plot.
#'
#' @seealso \code{\link{Incremental}}
#'
#' @export
plot.incremental <- function(x,
                             quantiles = c(0.05, 0.95),
                             ylab = "Performance",
                             pch = 18,
                             col = c("red", "grey"), col.axis = NULL,
                             cex = 1, cex.lab = 1.5,
                             xcex.axis = 1, ycex.axis = 1,
                             xlas = 2, ylas = 1,
                             sfrac = 0.005,
                             ylim = NULL, bty = "o",
                             xgrid = FALSE, ygrid = FALSE,
                             output_data = FALSE,
                             ...) {
  # Checking plotrix package is installed
  CheckPackageInstalled("plotrix")

  # Checking the inputs
  quantiles <- sort(quantiles)

  # Re-formatting the inputs
  if (is.null(col.axis)) {
    col.axis <- col
  }

  if ("concordance" %in% names(x)) {
    if ("lower" %in% names(x)) {
      xfull <- x$concordance
      xlower <- x$lower
      xupper <- x$upper
    } else {
      xfull <- sapply(x$concordance, stats::median, na.rm = TRUE)
      xlower <- sapply(x$concordance, stats::quantile, probs = quantiles[1], na.rm = TRUE)
      xupper <- sapply(x$concordance, stats::quantile, probs = quantiles[2], na.rm = TRUE)
    }
  }

  if ("AUC" %in% names(x)) {
    xfull <- sapply(x$AUC, stats::median, na.rm = TRUE)
    xlower <- sapply(x$AUC, stats::quantile, probs = quantiles[1], na.rm = TRUE)
    xupper <- sapply(x$AUC, stats::quantile, probs = quantiles[2], na.rm = TRUE)
  }

  if ("Q_squared" %in% names(x)) {
    xfull <- sapply(x$Q_squared, stats::median, na.rm = TRUE)
    xlower <- sapply(x$Q_squared, stats::quantile, probs = quantiles[1], na.rm = TRUE)
    xupper <- sapply(x$Q_squared, stats::quantile, probs = quantiles[2], na.rm = TRUE)
  }
  xseq <- 1:length(xfull)

  # Re-formatting label colours
  if (length(col.axis) < length(xfull)) {
    col.axis <- rep(col.axis, length(xfull))
  }

  # Defining the plot range
  if (is.null(ylim)) {
    ylim <- range(c(xlower, xupper, xfull))
  }

  # Defining horizontal grid
  hseq <- NULL
  if (ygrid) {
    hseq <- grDevices::axisTicks(ylim, log = FALSE)
  }

  # Defining vertical grid
  vseq <- NULL
  if (xgrid) {
    vseq <- xseq
  }

  # Defining colours by stable selection status
  if ("stable" %in% names(x)) {
    mycolours <- col[abs(x$stable - 2)]
  } else {
    mycolours <- col[1]
  }
  if ("stable" %in% names(x)) {
    mycolours_axis <- col.axis[abs(x$stable - 2)]
  } else {
    mycolours_axis <- col.axis[1]
  }

  # Creating the plot
  plot(NULL,
    cex.lab = cex.lab,
    bty = bty,
    xlim = range(xseq), ylim = ylim,
    panel.first = c(
      graphics::abline(h = hseq, col = "grey", lty = 3),
      graphics::abline(v = vseq, col = "grey", lty = 3)
    ),
    xlab = "", xaxt = "n", yaxt = "n", ylab = ylab,
    ...
  )
  graphics::axis(
    side = 2, at = grDevices::axisTicks(ylim, log = FALSE),
    las = ylas, cex.axis = ycex.axis
  )
  plotrix::plotCI(
    x = xseq, y = xfull, li = xlower, ui = xupper,
    pch = pch, cex = cex,
    sfrac = sfrac,
    col = mycolours, add = TRUE
  )
  graphics::axis(side = 1, at = xseq, labels = NA)
  for (k in 1:length(xseq)) {
    graphics::axis(
      side = 1, at = xseq[k],
      labels = ifelse(k == 1, yes = x$names[k], no = paste0("+ ", x$names[k])),
      col.axis = mycolours_axis[k],
      las = xlas, cex.axis = xcex.axis
    )
  }

  if (output_data) {
    mat <- rbind(xfull, xlower, xupper)
    colnames(mat) <- x$names
    return(mat)
  }
}

#' @rdname plot.incremental
#' @export
IncrementalPlot <- plot.incremental

#' @rdname plot.incremental
#' @export
PlotIncremental <- plot.incremental

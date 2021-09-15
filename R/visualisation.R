#' Calibration plot
#'
#' Creates a plot showing the stability score as a function of the parameter
#' controlling the level of sparsity in the underlying feature selection
#' algorithm and/or the threshold in selection proportions.
#'
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}} or \code{\link{Clustering}}.
#' @param metric parameter to visualise. Possible values are "lambda" (parameter
#'   controlling the level of sparsity in underlying algorithm), "pi" (threshold
#'   in selection proportion) or "both".
#' @param clustering logical indicating whether calibration for clustering
#'   (\code{clustering=TRUE}) or variable selection (\code{clustering=FALSE})
#'   should be represented. This argument is only used if \code{stability} is the
#'   output from \code{\link{Clustering}}.
#' @param block_id ID of the block to visualise. Only used for multi-block
#'   stability selection graphical models. With block_id=NULL, all blocks are
#'   represented in separate panels.
#' @param lines logical indicating if the points should be linked by lines.
#' @param colours vector of colours used for the heatmap. By default a gradient
#'   of colours ranging from ivory to dark red is used. Only used with
#'   metric="both".
#' @param legend logical indicating if the colour bar should be included. Only
#'   used with metric="both".
#' @param legend_length length of the colour bar. Only used with metric="both".
#' @param legend_range range of the colour bar. Only used with metric="both".
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @param zlab label of the z-axis.
#' @param filename file path to saved figure. With filename=NULL, the plot is
#'   not saved.
#' @param fileformat format of the saved figure. Possible values are "pdf" or
#'   "png". Only used if argument "filename" is not NULL.
#' @param res resolution of the png figure (see \code{\link{png}} from the
#'   grDevices package). Only used if argument "filename" is not NULL and
#'   fileformat="png".
#' @param width width of the saved figure. Only used if argument "filename" is
#'   not NULL.
#' @param height height of the saved figure. Only used if argument "filename" is
#'   not NULL.
#' @param units units of width and height. Possible values are "px", "in", "cm"
#'   and "mm" (see \code{\link{png}} from the grDevices package). Only used if
#'   argument "filename" is not NULL and fileformat="png".
#' @param ... additional arguments to be passed to \code{\link{Heatmap}} or
#'   \code{\link{axis}} from the graphics package.
#'
#' @return a calibration plot.
#'
#' @details When selecting a single parameter, each point represents the best
#'   (maximum) stability score across all visited values of the other parameter.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20, nu_within = 0.1)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#'
#' # Calibration plots
#' par(mar = c(7, 5, 7, 6))
#' CalibrationPlot(stab)
#' par(mar = c(7, 5, 7, 1))
#' CalibrationPlot(stab, metric = "lambda")
#' CalibrationPlot(stab, metric = "pi")
#'
#' # User-defined colours (heatmap)
#' par(mar = c(7, 5, 7, 6))
#' CalibrationPlot(stab, colours = c("lightgrey", "blue"))
#' CalibrationPlot(stab, colours = c("lightgrey", "blue", "black"))
#' CalibrationPlot(stab,
#'   colours = c("lightgrey", "blue", "black"),
#'   legend_length = 25, legend_range = c(0, 4000)
#' )
#' }
#'
#' @export
CalibrationPlot <- function(stability, clustering = FALSE, metric = "both", block_id = NULL,
                            lines = FALSE, colours = c("ivory", "navajowhite", "tomato", "darkred"),
                            legend = TRUE, legend_length = 15, legend_range = NULL,
                            xlab = expression(lambda), ylab = expression(pi), zlab = expression(italic(q)),
                            filename = NULL, fileformat = "pdf", res = 500,
                            width = 7, height = 7, units = "in", ...) {
  # Extracting the number of blocks/components
  if ((stability$methods$type == "graphical_model") & (is.null(block_id))) {
    bigblocks <- BlockMatrix(stability$params$pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- unique(as.vector(bigblocks_vect))
    names(N_blocks) <- blocks
    nblocks <- max(blocks)
    block_id <- 1:nblocks
  } else {
    block_id <- 1
  }
  nblocks <- length(block_id)

  # Saving as PDF
  if (!is.null(filename)) {
    if (fileformat == "pdf") {
      grDevices::pdf(filename, width = width, height = height)
    } else {
      grDevices::png(filename, width = width, height = height, res = res, units = units)
    }
  }

  if (metric == "both") {
    for (b in block_id) {
      # Extracting the stability scores
      if (clustering) {
        mat <- stability$Sc_2d
        if (length(unique(stability$Lambda)) > 1) {
          # Identifying best number of contributing variables
          lambda_hat <- stability$Lambda[which.max(stability$S), 1]
          ids <- which(as.character(stability$Lambda) == lambda_hat)
        } else {
          ids <- 1:nrow(stability$Sc)
        }
        mat <- mat[ids, ]
      } else {
        if (length(stability$params$pk) == 1) {
          mat <- stability$S_2d
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        } else {
          mat <- stability$S_2d[, , b]
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        }
      }
      mat <- mat[, , drop = FALSE]
      colnames(mat) <- stability$params$pi_list
      if (stability$methods$type == "clustering") {
        if (length(unique(stability$Lambda[, b])) > 1) {
          rownames(mat) <- paste0(stability$nc[, b], " - ", stability$Lambda[, b])[ids]
        } else {
          rownames(mat) <- (stability$nc[, b])[ids]
        }
      } else {
        if (grepl("penalised", tolower(stability$methods$implementation))) {
          rownames(mat) <- formatC(stability$Lambda[, b], format = "e", digits = 2)[ids]
        } else {
          rownames(mat) <- (stability$Lambda[, b])[ids]
        }
      }

      # Extracting corresponding numbers of selected variables (q)
      Q <- stability$Q[, b]
      Q <- Q[ids]

      # Heatmap representation
      Heatmap(t(mat[nrow(mat):1, ncol(mat):1]),
        colours = colours, axes = FALSE,
        legend = legend, legend_length = legend_length, legend_range = legend_range
      )

      # Identifying best pair of parameters
      withr::local_par(list(xpd = FALSE))
      if (stability$methods$type == "clustering") {
        if (clustering) {
          graphics::abline(v = nrow(mat) - which(stability$nc[ids, b] == Argmax(stability, clustering = clustering)[b, 1]) + 0.5, lty = 3)
        } else {
          tmp <- paste0(stability$nc[, b], " - ", stability$Lambda[, b])[ArgmaxId(stability, clustering = clustering)[1, 1]]
          graphics::abline(v = nrow(mat) - which(rownames(mat) == tmp) + 0.5, lty = 3)
        }
      } else {
        graphics::abline(v = nrow(mat) - which(stability$Lambda[ids, b] == Argmax(stability, clustering = clustering)[b, 1]) + 0.5, lty = 3)
      }
      graphics::abline(h = which.min(abs(as.numeric(colnames(mat)) - Argmax(stability, clustering = clustering)[b, 2])) - 0.5, lty = 3)

      # Including axes
      graphics::axis(
        side = 2, at = (1:ncol(mat)) - 0.5, las = 2,
        labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2), ...
      )
      if (grepl("penalised", tolower(stability$methods$implementation))) {
        graphics::axis(
          side = 3, at = (1:nrow(mat)) - 0.5, las = 2,
          labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0)), ...
        )
        graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)), ...)
      } else {
        graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)), ...)
      }

      # Including axis labels
      graphics::mtext(text = ylab, side = 2, line = 3.5, cex = 1.5)
      if (grepl("penalised", tolower(stability$methods$implementation))) {
        graphics::mtext(text = xlab, side = 1, line = 5.2, cex = 1.5)
        graphics::mtext(text = zlab, side = 3, line = 3.5, cex = 1.5)
      } else {
        graphics::mtext(text = xlab, side = 1, line = 3.5, cex = 1.5)
      }
    }
  } else {
    if (metric == "lambda") {
      for (b in block_id) {
        # Extracting the stability scores
        if (length(stability$params$pk) == 1) {
          mat <- stability$S_2d
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        } else {
          mat <- stability$S_2d[, , b]
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        }

        # Extracting the best stability score (with optimal pi) for each lambda value
        vect <- apply(mat, 1, max, na.rm = TRUE)

        # Extracting corresponding numbers of selected variables (q)
        Q <- stability$Q[, b, drop = FALSE]
        Q <- Q[ids]

        # Extracting corresponding lambda values
        Lambda <- stability$Lambda[ids, b, drop = FALSE]

        # Re-ordering by decreasing lambda
        ids <- sort.list(Lambda, decreasing = TRUE)
        Lambda <- Lambda[ids]
        Q <- Q[ids]
        vect <- vect[ids]

        # Using input ylab if not as default
        ylab <- ifelse(as.character(ylab) == as.character(expression(pi)), yes = "Stability Score", no = ylab)

        # Making plot
        cex_points <- 0.7
        plot(Lambda, vect,
          pch = 19, col = "navy", cex = cex_points,
          xlab = "", ylab = ylab, cex.lab = 1.5, xaxt = "n"
        )
        graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
        graphics::abline(h = max(vect), lty = 2, col = "red")
        graphics::abline(v = which.max(vect), lty = 2, col = "red")
        graphics::points(Lambda, vect, pch = 19, col = "navy", cex = cex_points)
        if (lines) {
          graphics::lines(Lambda, vect, col = "navy")
        }

        # Adding x-axis and z-axis and their labels
        xseq <- seq(1, length(Lambda), length.out = 5)
        graphics::abline(v = Lambda[xseq], lty = 3, col = "grey")
        graphics::axis(side = 1, at = Lambda[xseq], labels = formatC(Lambda[xseq], format = "e", digits = 2), las = 2)
        graphics::axis(
          side = 3, at = Lambda[xseq], las = 2,
          labels = formatC(Q[xseq], format = "f", big.mark = ",", digits = 0)
        )
        graphics::mtext(text = xlab, side = 1, line = 5.2, cex = 1.5)
        graphics::mtext(text = zlab, side = 3, line = 3.5, cex = 1.5)
      }
    }

    if (metric == "pi") {
      for (b in block_id) {
        # Extracting the stability scores
        if (length(stability$params$pk) == 1) {
          mat <- stability$S_2d
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        } else {
          mat <- stability$S_2d[, , b]
          ids <- which(apply(mat, 1, FUN = function(x) {
            any(!is.na(x))
          }))
          mat <- mat[ids, , drop = FALSE]
        }

        # Extracting the best stability score (with optimal lambda) for each pi value
        vect <- apply(mat, 2, max, na.rm = TRUE)

        # Using input ylab if not as default
        ylab <- ifelse(as.character(ylab) == as.character(expression(pi)), yes = "Stability Score", no = ylab)

        # Using input xlab if not as default
        xlab <- ifelse(as.character(xlab) == as.character(expression(lambda)), yes = expression(pi), no = ylab)

        # Making plot
        cex_points <- 0.7
        plot(1:length(vect), vect,
          pch = 19, col = "navy", cex = cex_points,
          xlab = "", ylab = ylab, cex.lab = 1.5, xaxt = "n"
        )
        xticks <- graphics::axTicks(side = 1)
        if (min(xticks) == 0) {
          xticks <- xticks + 1
        }
        graphics::abline(v = xticks, lty = 3, col = "grey")
        graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
        graphics::abline(h = max(vect), lty = 2, col = "red")
        graphics::abline(v = which.max(vect), lty = 2, col = "red")
        graphics::points(1:length(vect), vect, pch = 19, col = "navy", cex = cex_points)
        if (lines) {
          graphics::lines(1:length(vect), vect, col = "navy")
        }

        # Adding x-axis and its labels
        graphics::axis(side = 1, at = xticks, labels = formatC(stability$params$pi_list[xticks], digits = 2), las = 2)
        graphics::mtext(text = xlab, side = 1, line = 5.2, cex = 1.5)
      }
    }
  }

  if (!is.null(filename)) {
    grDevices::dev.off()
  }
}


#' Heatmap visualisation
#'
#' Produces a heatmap for visualisation of matrix entries.
#'
#' @param mat data matrix.
#' @param colours vector of colours used for the heatmap. By default a gradient
#'   of colours ranging from ivory to dark red is used.
#' @param resolution number of different colours to use.
#' @param axes logical indicating if the row and column names of \code{mat}
#'   should be displayed.
#' @param legend logical indicating if the colour bar should be included.
#' @param legend_length length of the colour bar.
#' @param legend_range range of the colour bar.
#'
#' @return A heatmap.
#'
#' @seealso \code{\link{CalibrationPlot}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' mat <- matrix(rnorm(200), ncol = 20)
#'
#' # Generating heatmaps
#' Heatmap(mat = mat)
#' Heatmap(mat = mat, colours = c("lightgrey", "blue", "black"), legend = FALSE)
#' }
#'
#' @export
Heatmap <- function(mat, colours = c("ivory", "navajowhite", "tomato", "darkred"),
                    resolution = 10000, axes = TRUE,
                    legend = TRUE, legend_length = NULL, legend_range = NULL) {
  # Transposing the input matrix so that rows are rows
  mat <- t(mat)

  # Defining the legend length
  if (is.null(legend_length)) {
    legend_length <- ncol(mat)
  }

  # Preparing colours
  colours <- grDevices::colorRampPalette(colours)(resolution)
  names(colours) <- 1:resolution

  # Re-formatting matrix
  mat <- mat[, ncol(mat):1]
  vect <- as.vector(mat)

  # Defining extreme values
  if (is.null(legend_range)) {
    # myrange <- c(min(vect, na.rm = TRUE), max(vect, na.rm = TRUE))
    myrange <- range(vect, na.rm = TRUE)
    myrange <- c(floor(myrange[1]), ceiling(myrange[2]))
  } else {
    myrange <- legend_range
  }

  # Getting corresponding colours
  mycol <- as.character(cut(vect, breaks = seq(myrange[1], myrange[2], length.out = resolution + 1), labels = 1:resolution, include.lowest = TRUE))
  mycol_mat <- matrix(mycol, ncol = ncol(mat))

  # Making heatmap
  plot(NA,
    xlim = c(0, nrow(mycol_mat)), ylim = c(0, ncol(mycol_mat)),
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
  )
  for (i in 0:(nrow(mycol_mat) - 1)) {
    for (j in 0:(ncol(mycol_mat) - 1)) {
      graphics::polygon(
        x = c(i, i + 1, i + 1, i), y = c(j, j, j + 1, j + 1),
        col = colours[mycol_mat[i + 1, j + 1]],
        border = colours[mycol_mat[i + 1, j + 1]]
      )
    }
  }
  if (axes) {
    if (!is.null(rownames(mat))) {
      graphics::axis(side = 1, at = 1:nrow(mat) - 0.5, labels = rownames(mat), las = 2)
    } else {
      graphics::axis(side = 1, at = c(1, nrow(mat)) - 0.5, labels = NA, las = 2)
    }
    if (!is.null(colnames(mat))) {
      graphics::axis(side = 2, at = 1:ncol(mat) - 0.5, labels = colnames(mat), las = 2)
    } else {
      graphics::axis(side = 2, at = c(1, ncol(mat)) - 0.5, labels = NA, las = 2)
    }
  }


  # Adding colour bar (legend)
  if (legend) {
    withr::local_par(list(xpd = TRUE))
    legend_width_factor <- 1.05
    mylegend_values <- grDevices::axisTicks(c(myrange[1], myrange[2]), log = FALSE)
    mylegend_ids <- as.numeric(as.character(cut(mylegend_values,
      breaks = seq(myrange[1], myrange[2], length.out = resolution + 1),
      labels = 1:resolution, include.lowest = TRUE
    )))
    ypos <- ncol(mat)
    xpos <- nrow(mat) * 1.05
    for (l in 1:length(colours)) {
      graphics::polygon(
        x = c(xpos, xpos * legend_width_factor, xpos * legend_width_factor, xpos),
        y = c(
          ypos - legend_length + legend_length * l / length(colours),
          ypos - legend_length + legend_length * l / length(colours),
          ypos - legend_length + legend_length * (l + 1) / length(colours),
          ypos - legend_length + legend_length * (l + 1) / length(colours)
        ),
        col = colours[l], border = colours[l]
      )
      if (l %in% mylegend_ids) {
        graphics::text(
          x = xpos * legend_width_factor, y = ypos - legend_length + legend_length * (l + 0.5) / length(colours),
          labels = paste0("- ", mylegend_values[which(mylegend_ids == l)]), adj = c(0, 0.5)
        )
      }
    }
    withr::local_par(list(xpd = FALSE)) # for legend
  }
}

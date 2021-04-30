#' Calibration plot
#'
#' Returns a heatmap or plot showing the stability score as a function of the
#' parameter controlling the sparsity and/or thereshold in selection proportion.
#'
#' @param stability output from \code{\link{VariableSelection}} or
#'   \code{\link{GraphicalModel}}.
#' @param metric parameter to visualise. Possible values are "lambda" (parameter
#'   controlling the level of sparsity in underlying algorithm), "pi" (threshold
#'   in selection proportion) or "both".
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
#' @param mar vector of margins (see \code{\link{par}} from the graphics
#'   package). With mar=NULL, margins are automatically defined.
#' @param mfrow vector defining the layout of the figure (see \code{\link{par}}
#'   from the graphics package). With mfrow=NULL, the generated figure has as
#'   many panels as there blocks in the model.
#' @param ... additional arguments to be passed to \code{\link{Heatmap}} or
#'   \code{\link{axis}} from the graphics package.
#'
#' @return a calibration plot.
#'
#' @details When selecting a single parameter, each point represents the best
#' (maximum) stability score across all visited values of the other parameter.
#'
#' @family calibration functions
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20, nu = 0.1)
#'
#' # Stability selection
#' stab <- GraphicalModel(data = simul$data)
#'
#' # Calibration plots
#' CalibrationPlot(stab)
#' CalibrationPlot(stab, metric = "lambda")
#' CalibrationPlot(stab, metric = "pi")
#'
#' # User-defined colours (heatmap)
#' CalibrationPlot(stab, colours = c("lightgrey", "blue"))
#' CalibrationPlot(stab, colours = c("lightgrey", "blue", "black"))
#' CalibrationPlot(stab,
#'   colours = c("lightgrey", "blue", "black"),
#'   legend_length = 25, legend_range = c(0, 4000)
#' )
#'
#' @export
CalibrationPlot <- function(stability, metric = "both", block_id = NULL,
                            lines = TRUE, colours = c("ivory", "navajowhite", "tomato", "darkred"),
                            legend = TRUE, legend_length = 15, legend_range = NULL,
                            xlab = expression(lambda), ylab = expression(pi), zlab = expression(italic(q)),
                            filename = NULL, fileformat = "pdf", res = 500,
                            width = 7, height = 7, units = "in", mar = NULL, mfrow = NULL, ...) {
  # Extracting the number of blocks
  if (is.null(block_id)) {
    bigblocks <- BlockMatrix(stability$params$pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- unique(as.vector(bigblocks_vect))
    names(N_blocks) <- blocks
    nblocks <- max(blocks)
    block_id <- 1:nblocks
  } else {
    nblocks <- 1
  }

  # Saving as PDF
  if (!is.null(filename)) {
    if (fileformat == "pdf") {
      grDevices::pdf(filename, width = width, height = height)
    } else {
      grDevices::png(filename, width = width, height = height, res = res, units = units)
    }
  }

  if (metric == "both") {
    if (is.null(mar)) {
      mar <- c(7, 5, 7, 7)
    }
    if (is.null(mfrow)) {
      mfrow <- c(1, nblocks)
    }
    graphics::par(mar = mar, mfrow = mfrow)

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
      mat <- mat[, , drop = FALSE]
      colnames(mat) <- stability$params$pi_list
      rownames(mat) <- formatC(stability$Lambda[, b], format = "e", digits = 2)[ids]

      # Extracting corresponding numbers of selected variables (q)
      Q <- stability$Q[, b]
      Q <- Q[ids]

      # Heatmap representation
      Heatmap(mat[nrow(mat):1, ncol(mat):1],
        colours = colours,
        legend = legend, legend_length = legend_length, legend_range = legend_range
      )

      # Identifying best pair of parameters
      graphics::abline(h = which.min(abs(as.numeric(colnames(mat)) - Argmax(stability)[b, 2])) - 0.5, lty = 3)
      graphics::abline(v = nrow(mat) - which(stability$Lambda[ids, b] == Argmax(stability)[b, 1]) + 0.5, lty = 3)

      # Including axes
      graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)), ...)
      graphics::axis(
        side = 2, at = (1:ncol(mat)) - 0.5, las = 2,
        labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2), ...
      )
      graphics::axis(
        side = 3, at = (1:nrow(mat)) - 0.5, las = 2,
        labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0)), ...
      )

      # Including axis labels
      graphics::mtext(text = xlab, side = 1, line = 5.2, cex = 1.5)
      graphics::mtext(text = ylab, side = 2, line = 3.5, cex = 1.5)
      graphics::mtext(text = zlab, side = 3, line = 3.5, cex = 1.5)
    }
  } else {
    if (metric == "lambda") {
      if (is.null(mar)) {
        mar <- c(7, 5, 7, 2)
      }
      if (is.null(mfrow)) {
        mfrow <- c(1, nblocks)
      }
      graphics::par(mar = mar, mfrow = mfrow)

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
      if (is.null(mar)) {
        mar <- c(7, 5, 2, 2)
      }
      if (is.null(mfrow)) {
        mfrow <- c(1, nblocks)
      }
      graphics::par(mar = mar, mfrow = mfrow)

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


#' Heatmap
#'
#' Produces a heatmap from a matrix.
#'
#' @param mat data matrix.
#' @param colours vector of colours used for the heatmap.
#' By default a gradient of colours ranging from ivory to dark red is used.
#' @param resolution number of different colours to use.
#' @param legend logical indicating if the colour bar should be included.
#' @param legend_length length of the colour bar.
#' @param legend_range range of the colour bar.
#'
#' @return a heatmap.
#'
#' @seealso \code{\link{CalibrationPlot}}
#'
#' @examples
#' # Data simulation
#' set.seed(1)
#' mat <- matrix(rnorm(200), ncol = 20)
#'
#' # Generating heatmaps
#' Heatmap(mat = mat)
#' Heatmap(mat = mat, colours = c("lightgrey", "blue", "black"), legend = FALSE)
#'
#' @export
Heatmap <- function(mat, colours = c("ivory", "navajowhite", "tomato", "darkred"),
                    resolution = 10000, legend = TRUE, legend_length = 15, legend_range = NULL) {
  # Preparing colours
  colours <- grDevices::colorRampPalette(colours)(resolution)
  names(colours) <- 1:resolution

  # Re-formatting matrix
  mat <- mat[, ncol(mat):1]
  vect <- as.vector(mat)

  # Defining extreme values
  if (is.null(legend_range)) {
    myrange <- c(min(vect, na.rm = TRUE), max(vect, na.rm = TRUE))
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

  # Adding colour bar (legend)
  if (legend) {
    graphics::par(xpd = TRUE)
    legend_width_factor <- 1.05
    myrange <- round(myrange)
    mylegend_values <- seq(myrange[1], myrange[2], length.out = 100)
    mylegend_values <- unique(round(mylegend_values, digits = -max(nchar(round(myrange))) + 2))
    if (is.null(legend_range)) {
      mylegend_values <- mylegend_values[mylegend_values >= min(myrange)]
    }
    if (length(mylegend_values) > legend_length) {
      mylegend_values <- unique(round(mylegend_values, digits = -max(nchar(round(myrange))) + 1))
    }
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
    graphics::par(xpd = FALSE)
  }
}

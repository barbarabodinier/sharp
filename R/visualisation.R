#' Calibration plot
#'
#' Creates a plot showing the stability score as a function of the parameter
#' controlling the level of sparsity in the underlying feature selection
#' algorithm and/or the threshold in selection proportions.
#'
#' @param stability output of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}} or \code{\link{BiSelection}}.
#' @param block_id ID of the block to visualise. Only used for multi-block
#'   stability selection graphical models. If \code{block_id=NULL}, all blocks
#'   are represented in separate panels.
#' @param heatmap logical indicating if results should be visualised as a
#'   heatmap. If \code{heatmap=TRUE}, the stability score is colour-coded and
#'   represented as a function of the parameter controlling the level of
#'   sparsity (x-axis) and threshold in selection proportions (y-axis). If
#'   \code{heatmap=FALSE}, the stability score (y-axis) is represented as a
#'   function of the parameter controlling the level of sparsity (x-axis). Only
#'   used if \code{stability} is of class \code{graphical_model} or
#'   \code{variable_selection}.
#' @param colours vector of colours.
#' @param pch type of point, as in \code{\link{points}}.
#' @param cex size of point.
#' @param xlim displayed range along the x-axis.
#' @param ylim displayed range along the y-axis.
#' @param bty character string indicating if the box around the plot should be
#'   drawn. Possible values include: \code{"o"} (default, the box is drawn), or
#'   \code{"n"} (no box).
#' @param lines logical indicating if the points should be linked by lines.
#' @param lty line type for the point-wise median curve, as in
#'   \code{\link{par}}.
#' @param lwd line width for the point-wise median curve, as in
#'   \code{\link{par}}.
#' @param show_argmax logical indicating if the calibrated parameter(s) should
#'   be indicated by lines.
#' @param show_pix logical indicating if the calibrated threshold in selection
#'   proportion in \code{X} should be written for each point. Only used if
#'   \code{stability} is of class \code{bi_selection}.
#' @param show_piy logical indicating if the calibrated threshold in selection
#'   proportion in \code{Y} should be written for each point. Only used if
#'   \code{stability} is of class \code{bi_selection} with sparsity on a
#'   multivariate outcome.
#' @param offset distance between the point and the text, as in
#'   \code{\link[graphics]{text}}.
#' @param legend logical indicating if the legend should be included. Only used
#'   if \code{heatmap=TRUE} or with multiple components if \code{stability} is
#'   of class \code{bi_selection}.
#' @param legend_length length of the colour bar. Only used if
#'   \code{heatmap=TRUE}.
#' @param legend_range range of the colour bar. Only used if
#'   \code{heatmap=TRUE}.
#' @param xlab label of the x-axis.
#' @param ylab label of the y-axis.
#' @param zlab label of the z-axis.
#' @param xlas orientation of labels on the x-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylas orientation of labels on the y-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param cex.lab size of label on the y-axis.
#' @param cex.axis size of labels along the axes.
#' @param xgrid logical indicating if a vertical grid should be drawn.
#' @param ygrid logical indicating if a horizontal grid should be drawn.
#' @param params vector of possible parameters if \code{stability} is of class
#'   \code{bi_selection}. The order of these parameters defines the order in
#'   which they are represented.
#'
#' @return a calibration plot.
#'
#' @details When selecting a single parameter, each point represents the best
#'   (maximum) stability score across all visited values of the other parameter.
#'
#' @family calibration functions
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}},
#'   \code{\link{BiSelection}}
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
#' # Calibration heatmap
#' par(mar = c(7, 5, 6, 6))
#' CalibrationPlot(stab)
#'
#' # Calibration curve
#' par(mar = c(7, 5, 6, 1))
#' CalibrationPlot(stab, heatmap = FALSE)
#'
#' # User-defined colours
#' par(mar = c(7, 5, 7, 6))
#' CalibrationPlot(stab,
#'   colours = c("ivory", "blue", "black"),
#'   legend_length = 31,
#'   legend_range = c(0, 2500)
#' )
#' }
#'
#' @export
CalibrationPlot <- function(stability, block_id = NULL,
                            heatmap = TRUE, colours = NULL,
                            pch = 19, cex = 0.7,
                            xlim = NULL, ylim = NULL, bty = "o",
                            lines = TRUE, lty = 3, lwd = 2,
                            show_argmax = TRUE,
                            show_pix = FALSE, show_piy = FALSE, offset = 0.3,
                            legend = TRUE, legend_length = 15, legend_range = NULL,
                            xlab = NULL, ylab = NULL, zlab = expression(italic(q)),
                            xlas = 2, ylas = 0, cex.lab = 1.5, cex.axis = 1,
                            xgrid = FALSE, ygrid = FALSE,
                            params = c("ny", "alphay", "nx", "alphax")) {
  # To deal with later: showing calibration of clustering or selection
  clustering <- FALSE

  if (class(stability) == "bi_selection") {
    # Extracting summary information
    x <- stability$summary_full

    # Checking input
    params <- unique(params)
    if (any(!c("ny", "alphay", "nx", "alphax") %in% params)) {
      params <- c("ny", "alphay", "nx", "alphax")
      warning("Invalid input for argument 'params'. Please provide a vector with all the following: 'ny', 'alphay', 'nx', 'alphax'.")
    }
    params <- params[params %in% c("ny", "alphay", "nx", "alphax")]

    # Identifying parameters
    params <- params[params %in% colnames(x)]
    print(params)

    # Defining default arguments
    if (is.null(ylab)) {
      ylab <- "Stability Score"
    }

    if (is.null(xlab)) {
      if (length(params) > 1) {
        xlab <- ""
      } else {
        xlab <- expression(n[X])
      }
    }

    if (is.null(colours)) {
      colours <- grDevices::colorRampPalette(c("navy", "darkred"))(nrow(stability$summary))
    } else {
      colours <- grDevices::colorRampPalette(colours)(nrow(stability$summary))
    }

    if (length(unique(x$comp)) == 1) {
      legend <- FALSE
    }

    if (is.null(xlim)) {
      xlim <- c(0.5, max(sapply(split(x, f = x$comp), nrow)) + 0.5)
    }

    if (is.null(ylim)) {
      ylim <- range(x$S)
      if (legend) {
        ylim[2] <- ylim[2] + diff(ylim) * 0.15
      }
    }

    # Drawing one set of points per component
    for (comp_id in unique(x$comp)) {
      tmp <- x[which(x$comp == comp_id), ]

      # Ensuring increasing ny
      tmp <- tmp[do.call(order, tmp[, params, drop = FALSE]), ]
      # if ("ny" %in% colnames(tmp)) {
      #   tmp=tmp[order(tmp$ny, tmp$nx), ]
      # }
      # tmp=tmp[order(lapply(params, FUN=function(param_id){with(tmp, eval(parse(text=param_id)))})),]
      # if ("alphax" %in% colnames(tmp))

      if (comp_id == min(x$comp)) {
        # Initialising the plot
        plot(NA,
          xlim = xlim, ylim = ylim, bty = bty,
          xlab = xlab, ylab = ylab, cex.lab = cex.lab,
          cex.axis = cex.axis,
          xaxt = "n", las = ylas
        )

        # Defining vertical grid
        if (xgrid) {
          withr::local_par(list(xpd = FALSE))
          graphics::abline(v = 1:nrow(tmp), lty = 3, col = "grey")
        }

        # Defining horizontal grid
        if (ygrid) {
          withr::local_par(list(xpd = FALSE))
          graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
        }

        # Adding x-axes
        for (param_id in 1:length(params)) {
          if (param_id == 1) {
            graphics::axis(
              side = 1, at = 1:nrow(tmp),
              labels = tmp[, rev(params)[param_id]],
              cex.axis = cex.axis, las = xlas
            )
          } else {
            ids <- c(1, which(diff(tmp[, rev(params)[param_id]]) != 0) + 1)
            ids <- c(ids - 0.5, nrow(tmp) + 0.5)
            graphics::axis(side = 1, at = ids, labels = NA, line = (param_id - 1) * 3)
            withr::local_par(list(xpd = FALSE))
            graphics::abline(v = ids, lty = 2)
            ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
            graphics::axis(
              side = 1, at = ids, labels = tmp[ids, rev(params)[param_id]],
              line = (param_id - 1) * 3, tick = FALSE, cex.axis = cex.axis, las = xlas
            )
          }
        }
        # graphics::axis(side = 1, at = 1:nrow(tmp), labels = tmp$nx, cex.axis = cex.axis, las = xlas)
        # if ("ny" %in% colnames(tmp)) {
        #   ids <- c(which(!duplicated(tmp$ny)) - 0.5, nrow(tmp) + 0.5)
        #   graphics::axis(side = 1, at = ids, labels = NA, line = 3)
        #   withr::local_par(list(xpd = FALSE))
        #   graphics::abline(v = ids, lty = 2)
        #   ids <- apply(rbind(ids[-1], ids[-length(ids)]), 2, mean)
        #   graphics::axis(side = 1, at = ids, labels = unique(tmp$ny), line = 3, tick = FALSE, cex.axis = cex.axis, las = xlas)
        # }

        # Adding x-labels
        if (length(params) > 1) {
          for (param_id in 1:length(params)) {
            if (rev(params)[param_id] == "nx") {
              graphics::mtext(text = expression(n[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
            }
            if (rev(params)[param_id] == "alphax") {
              graphics::mtext(text = expression(alpha[X]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
            }
            if (rev(params)[param_id] == "ny") {
              graphics::mtext(text = expression(n[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
            }
            if (rev(params)[param_id] == "alphay") {
              graphics::mtext(text = expression(alpha[Y]), side = 1, at = 0.5 - nrow(tmp) * 0.1, line = (param_id - 1) * 3 + 1, cex = cex.lab)
            }
          }
          # graphics::mtext(text = expression(n[X]), side = 1, at = -nrow(tmp) * 0.06, line = 1, cex = cex.lab)
          # graphics::mtext(text = expression(n[Y]), side = 1, at = -nrow(tmp) * 0.06, line = 4, cex = cex.lab)
        }
      }

      # Adding calibrated lines
      if (show_argmax) {
        withr::local_par(list(xpd = FALSE))
        graphics::abline(v = which.max(tmp$S), lty = 3, col = colours[comp_id])
      }

      # Adding lines
      if (lines) {
        # if ("ny" %in% colnames(tmp)) {
        #   for (y_value in unique(tmp$ny)) {
        #     graphics::lines(which(tmp$ny == y_value),
        #                     tmp[which(tmp$ny == y_value), "S"],
        #                     col = colours[comp_id],
        #                     lty = lty, lwd = lwd
        #     )
        #   }
        # } else {
        graphics::lines(1:nrow(tmp),
          tmp$S,
          col = colours[comp_id],
          lty = lty, lwd = lwd
        )
        # }
      }

      # Adding data points
      graphics::points(tmp$S,
        pch = pch,
        col = colours[comp_id],
        cex = cex
      )

      # Adding pi values
      if ((show_pix) & (!show_piy)) {
        graphics::text(tmp$S,
          labels = tmp$pix,
          col = colours[comp_id],
          cex = cex, pos = 3,
          offset = offset
        )
      }

      if ((!show_pix) & (show_piy)) {
        graphics::text(tmp$S,
          labels = tmp$piy,
          col = colours[comp_id],
          cex = cex, pos = 3,
          offset = offset
        )
      }

      if ((show_pix) & (show_piy)) {
        for (k in 1:nrow(tmp)) {
          graphics::text(k, tmp[k, "S"],
            labels = eval(parse(text = paste0("expression(pi[x]*' = ", tmp[k, "pix"], " ; '*pi[y]*' = ", tmp[k, "piy"], "')"))),
            col = colours[comp_id],
            cex = cex, pos = 3,
            offset = offset
          )
        }
      }
    }

    # Adding legend
    if (legend) {
      graphics::legend("top",
        col = colours, lty = lty, pch = pch, lwd = lwd,
        legend = paste0("Component ", unique(x$comp)),
        horiz = TRUE, bg = "white"
      )
    }
  } else {
    # Defining default arguments
    if (heatmap) {
      metric <- "both"
      if (is.null(colours)) {
        colours <- c("ivory", "navajowhite", "tomato", "darkred")
      }
      if (is.null(ylab)) {
        ylab <- expression(pi)
      }
    } else {
      metric <- "lambda"
      if (is.null(colours)) {
        colours <- "navy"
      }
      if (is.null(ylab)) {
        ylab <- "Stability Score"
      }
    }
    if (is.null(xlab)) {
      xlab <- expression(lambda)
    }

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
          colours = colours, bty = bty, axes = FALSE,
          legend = legend, legend_length = legend_length, legend_range = legend_range
        )

        # Adding calibrated lines
        if (show_argmax) {
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
        }

        # Including axes
        graphics::axis(
          side = 2, at = (1:ncol(mat)) - 0.5, las = 2, cex.axis = cex.axis,
          labels = formatC(as.numeric(colnames(mat)), format = "f", digits = 2)
        )
        if (grepl("penalised", tolower(stability$methods$implementation))) {
          graphics::axis(
            side = 3, at = (1:nrow(mat)) - 0.5, las = 2, cex.axis = cex.axis,
            labels = rev(formatC(Q, format = "f", big.mark = ",", digits = 0))
          )
          graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)), cex.axis = cex.axis)
        } else {
          graphics::axis(side = 1, at = (1:nrow(mat)) - 0.5, las = 2, labels = rev(rownames(mat)), cex.axis = cex.axis)
        }

        # Including axis labels
        graphics::mtext(text = ylab, side = 2, line = 3.5, cex = cex.lab)
        if (grepl("penalised", tolower(stability$methods$implementation))) {
          graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
          graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
        } else {
          graphics::mtext(text = xlab, side = 1, line = 3.5, cex = cex.lab)
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

          if (is.null(xlim)) {
            xlim <- range(Lambda, na.rm = TRUE)
          }

          if (is.null(ylim)) {
            ylim <- range(vect)
          }

          # Initialising the plot
          plot(NA,
            xlim = xlim, ylim = ylim, bty = bty,
            xlab = "", ylab = ylab, cex.lab = cex.lab,
            cex.axis = cex.axis,
            xaxt = "n", las = ylas
          )

          # Defining horizontal grid
          if (ygrid) {
            withr::local_par(list(xpd = FALSE))
            graphics::abline(h = graphics::axTicks(side = 2), lty = 3, col = "grey")
          }

          # Adding calibrated lines
          if (show_argmax) {
            withr::local_par(list(xpd = FALSE))
            graphics::abline(h = max(vect), lty = 3, col = colours[1])
            graphics::abline(v = Lambda[which.max(vect)], lty = 3, col = colours[1])
          }

          # Adding lines
          if (lines) {
            graphics::lines(Lambda, vect, col = colours[1], lty = lty, lwd = lwd)
          }

          # Adding data points
          graphics::points(Lambda, vect, pch = pch, col = colours[1], cex = cex)

          # Adding x-axis and z-axis and their labels
          lseq <- grDevices::axisTicks(range(Lambda, na.rm = TRUE), log = FALSE)
          xseq <- 1
          for (i in 1:length(lseq)) {
            xseq <- c(xseq, which.min(abs(Lambda - lseq[i])))
          }
          xseq <- c(xseq, length(Lambda))
          xseq <- unique(xseq)
          if (xgrid) {
            withr::local_par(list(xpd = FALSE))
            graphics::abline(v = Lambda[xseq], lty = 3, col = "grey")
          }
          graphics::axis(side = 1, at = Lambda[xseq], labels = formatC(Lambda[xseq], format = "e", digits = 2), las = xlas, cex.axis = cex.axis)
          graphics::axis(
            side = 3, at = Lambda[xseq], las = xlas,
            labels = formatC(Q[xseq], format = "f", big.mark = ",", digits = 0), cex.axis = cex.axis
          )
          graphics::mtext(text = xlab, side = 1, line = 5.2, cex = cex.lab)
          graphics::mtext(text = zlab, side = 3, line = 3.5, cex = cex.lab)
        }
      }
    }
  }
}


#' Heatmap visualisation
#'
#' Produces a heatmap for visualisation of matrix entries.
#'
#' @inheritParams CalibrationPlot
#' @param mat data matrix.
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
#' par(mar = c(1, 1, 1, 5))
#' Heatmap(mat = mat)
#' Heatmap(
#'   mat = mat,
#'   colours = c("lightgrey", "blue", "black"),
#'   legend = FALSE
#' )
#' }
#'
#' @export
Heatmap <- function(mat, colours = c("ivory", "navajowhite", "tomato", "darkred"),
                    resolution = 10000, bty = "o", axes = TRUE, cex.axis = 1,
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
  withr::local_par(xaxs = "i", yaxs = "i")
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
      graphics::axis(side = 1, at = 1:nrow(mat) - 0.5, labels = rownames(mat), las = 2, cex.axis = cex.axis)
    }
    if (!is.null(colnames(mat))) {
      graphics::axis(side = 2, at = 1:ncol(mat) - 0.5, labels = colnames(mat), las = 2, cex.axis = cex.axis)
    }
  }
  if (bty == "o") {
    graphics::box()
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

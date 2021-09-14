#' @export
print.variable_selection <- function(x, ...) {
  cat(paste0("Stability selection using function ", x$methods$implementation, " with ", x$methods$family, " family."))
  cat("\n")
  if (x$methods$resampling == "subsampling") {
    cat(paste0("The model was run using ", x$params$K, " subsamples of ", x$params$tau * 100, "% of the observations."))
  } else {
    cat(paste0("The model was run using ", x$params$K, " bootstrap samples."))
  }
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
      round(sum(object$selectedX[k, ]))
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
        round(sum(object$selectedY[k, ]))
      ))
      cat("\n")
    }
  }
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

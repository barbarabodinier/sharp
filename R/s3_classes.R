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
print.simulation_graphical_model <- function(x, ...) {
  cat(paste0("Multivariate Normal data with underlying structure of a graphical model."))
  cat("\n")
  cat("\n")
  cat(paste0("Number of observations: ", nrow(x$data)))
  cat("\n")
  cat(paste0("Number of variables (nodes): ", ncol(x$data)))
  cat("\n")
  cat(paste0("Number of edges: ", sum(x$theta == 1) / 2))
  cat("\n")
}


#' @export
print.simulation_clustering <- function(x, ...) {
  cat(paste0("Multivariate Normal data with underlying clusters of participants along (a subset of) variables."))
  cat("\n")
  cat("\n")
  cat(paste0("Number of observations: ", nrow(x$data)))
  cat("\n")
  cat(paste0("Number of clusters: ", max(x$theta)))
  for (k in 1:max(x$theta)) {
    cat("\n")
    cat(paste0("- Cluster ", k, " (N=", sum(x$theta == k), " observations)"))
  }
  cat("\n")
  cat("\n")
  cat(paste0("Number of variables: ", ncol(x$data)))
  cat("\n")
  cat(paste0("Number of variables contributing to the clustering: ", sum(x$theta_xc)))
  cat("\n")
}


#' @export
print.simulation_components <- function(x, ...) {
  cat(paste0("Multivariate Normal data with independent groups of variables."))
  cat("\n")
  cat("\n")
  cat(paste0("Number of observations: ", nrow(x$data)))
  cat("\n")
  cat("\n")
  cat(paste0("Number of variables: ", ncol(x$data)))
  cat("\n")
  cat(paste0("Number of independent groups of variables: ", max(x$membership)))
  for (k in 1:max(x$membership)) {
    cat("\n")
    cat(paste0("- Group ", k, " (N=", sum(x$membership == k), " variables)"))
  }
  cat("\n")
}


#' @export
print.simulation_regression <- function(x, ...) {
  cat(paste0("Multivariate Normal data with predictors and outcome(s)."))
  cat("\n")
  cat("\n")
  cat(paste0("Number of observations: ", nrow(x$xdata)))
  cat("\n")
  cat(paste0("Number of outcome variable(s): ", ncol(x$ydata)))
  cat("\n")
  cat(paste0("Number of predictor variables: ", ncol(x$xdata)))
  cat("\n")
  cat(paste0(
    "Number of predictor variables contributing to the outcome(s): ",
    sum(apply(x$beta, 1, sum) != 0)
  ))
  cat("\n")
}


#' @export
plot.simulation_graphical_model <- function(x, ...) {
  igraph::plot.igraph(Graph(x$theta), ...)
}


#' @export
plot.simulation_clustering <- function(x, ...) {
  # Visualisation of Euclidian distances along the contributing variable
  Heatmap(
    mat = as.matrix(stats::dist(x$data[, which(x$theta_xc == 1), drop = FALSE])),
    colours = c("navy", "white", "red")
  )
  graphics::title("Distances across variables contributing to clustering")
}


#' @export
plot.simulation_components <- function(x, ...) {
  Heatmap(
    mat = stats::cor(x$data),
    colours = c("navy", "white", "red"),
    legend_range = c(-1, 1)
  )
  graphics::title("Pearson's correlations")
}


#' @export
plot.simulation_regression <- function(x, ...) {
  plot(Graph(x$adjacency,
    satellites = TRUE,
    node_colour = c(
      rep("red", ncol(x$ydata)),
      rep("orange", ncol(x$zdata)),
      rep("skyblue", ncol(x$xdata))
    )
  ))
}

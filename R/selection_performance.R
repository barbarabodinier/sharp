#' Selection performance
#'
#' Computes different metrics of selection performance by comparing the set of
#' selected features to the set of true predictors/edges. This function can only
#' be used in simulation studies (i.e. when the true model is known).
#'
#' @inheritParams GraphicalModel
#' @param theta output from \code{\link{VariableSelection}},
#'   \code{\link{BiSelection}}, or \code{\link{GraphicalModel}}. Alternatively,
#'   it can be a binary matrix of selected variables (in variable selection) or
#'   a binary adjacency matrix (in graphical modelling)
#' @param theta_star output from \code{\link{SimulateRegression}},
#'   \code{\link{SimulateComponents}}, or \code{\link{SimulateGraphical}}.
#'   Alternatively, it can be a binary matrix of true predictors (in variable
#'   selection) or the true binary adjacency matrix (in graphical modelling).
#' @param cor optional correlation matrix. Only used in graphical modelling.
#' @param thr optional threshold in correlation. Only used in graphical
#'   modelling and when argument "cor" is not NULL.
#'
#' @return A matrix of selection metrics including:
#'
#'   \item{TP}{number of True Positives (TP)} \item{FN}{number of False
#'   Negatives (TN)} \item{FP}{number of False Positives (FP)} \item{TN}{number
#'   of True Negatives (TN)} \item{sensitivity}{sensitivity, i.e. TP/(TP+FN)}
#'   \item{specificity}{specificity, i.e. TN/(TN+FP)} \item{accuracy}{accuracy,
#'   i.e. (TP+TN)/(TP+TN+FP+FN)} \item{precision}{precision (p), i.e.
#'   TP/(TP+FP)} \item{recall}{recall (r), i.e. TP/(TP+FN)}
#'   \item{F1_score}{F1-score, i.e. 2*p*r/(p+r)}
#'
#'   If argument "cor" is provided, the number of False Positives among
#'   correlated (FP_c) and uncorrelated (FP_i) pairs, defined as having
#'   correlations (provided in "cor") above or below the threshold "thr", are
#'   also reported.
#'
#'   Block-specific performances are reported if "pk" is not NULL. In this case,
#'   the first row of the matrix corresponds to the overall performances, and
#'   subsequent rows correspond to each of the blocks. The order of the blocks
#'   is defined as in \code{\link{BlockStructure}}.
#'
#' @family functions for model performance
#'
#' @examples
#' \dontrun{
#'
#' # Variable selection model
#' set.seed(1)
#' simul <- SimulateRegression(pk = 30)
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' perf <- SelectionPerformance(theta = stab, theta_star = simul)
#' perf <- SelectionPerformance(
#'   theta = SelectedVariables(stab),
#'   theta_star = simul$theta
#' ) # alternative formulation
#'
#' # Single-block graphical model
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 30)
#' stab <- GraphicalModel(xdata = simul$data)
#' perf <- SelectionPerformance(theta = stab, theta_star = simul)
#' perf <- SelectionPerformance(
#'   theta = stab, theta_star = simul,
#'   cor = cor(simul$data), thr = 0.5
#' )
#' perf <- SelectionPerformance(
#'   theta = Adjacency(stab),
#'   theta_star = simul$theta
#' ) # alternative formulation
#'
#' # Multi-block graphical model
#' set.seed(1)
#' simul <- SimulateGraphical(pk = c(10, 10))
#' stab <- GraphicalModel(xdata = simul$data, pk = c(10, 10), lambda_other_blocks = rep(0, 3))
#' perf <- SelectionPerformance(theta = stab, theta_star = simul, pk = c(10, 10))
#' perf <- SelectionPerformance(
#'   theta = stab, theta_star = simul, pk = c(10, 10),
#'   cor = cor(simul$data), thr = 0.5
#' )
#' perf <- SelectionPerformance(
#'   theta = Adjacency(stab),
#'   theta_star = simul$theta,
#'   pk = c(10, 10)
#' ) # alternative formulation
#'
#' # Sparse PLS model
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#' stab <- BiSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   LambdaY = 1:(ncol(y) - 1),
#'   implementation = SparsePLS,
#'   n_cat = 2
#' )
#' perf <- SelectionPerformance(theta = stab, theta_star = simul)
#' perf <- SelectionPerformance(
#'   theta = stab$selected,
#'   theta_star = simul$theta
#' ) # alternative formulation
#' }
#'
#' @export
SelectionPerformance <- function(theta, theta_star, pk = NULL, cor = NULL, thr = 0.5) {
  # Re-formatting input theta
  if (any(class(theta) %in% c("variable_selection", "graphical_model", "bi_selection"))) {
    if (class(theta) == "graphical_model") {
      theta <- Adjacency(theta)
    } else {
      if (class(theta) == "variable_selection") {
        theta <- SelectedVariables(theta)
        theta <- as.vector(theta)
      } else {
        if ("selected" %in% names(theta)) {
          theta <- theta$selected # PLS
        } else {
          theta <- theta$selectedX # PCA
        }
      }
    }
  }

  # Re-formatting input theta_star
  if (any(class(theta_star) %in% c("simulation_regression", "simulation_graphical_model", "simulation_components"))) {
    theta_star <- theta_star$theta
  }
  if (is.vector(theta)) {
    theta_star <- as.vector(theta_star)
  } else {
    if (ncol(theta) != ncol(theta_star)) {
      theta_star <- theta_star[, 1:ncol(theta)]
    }
  }

  # Storing similarities/differences between estimated and true sets
  Asum <- theta + 2 * theta_star

  # Extracting block-specific performances
  if (is.null(pk)) {
    if (is.vector(Asum)) {
      out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
    } else {
      if (isSymmetric(theta_star)) {
        out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
      } else {
        out <- NULL
        for (k in 1:ncol(Asum)) {
          out <- rbind(out, SelectionPerformanceSingle(Asum[, k], cor = cor, thr = thr))
        }
        rownames(out) <- colnames(Asum)
      }
    }
  } else {
    Asum_vect <- Asum[upper.tri(Asum)]
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    if (!is.null(cor)) {
      cor_vect <- cor[upper.tri(cor)]
    } else {
      cor_vect <- NULL
    }

    out <- SelectionPerformanceSingle(Asum, cor = cor, thr = thr)
    for (k in sort(unique(bigblocks_vect))) {
      tmp <- SelectionPerformanceSingle(Asum_vect[bigblocks_vect == k],
        cor = cor_vect[bigblocks_vect == k], thr = thr
      )
      out <- rbind(out, tmp)
    }
  }

  return(out)
}


#' Edge-wise comparison of two graphs
#'
#' Generates an \code{\link[igraph:igraph-package]{igraph}} object representing
#' the common and graph-specific edges.
#'
#' @inheritParams Graph
#' @param graph1 first graph. Possible inputs are: adjacency matrix, or
#'   \code{\link[igraph:igraph-package]{igraph}} object, or output of
#'   \code{\link{GraphicalModel}}, \code{\link{VariableSelection}},
#'   \code{\link{BiSelection}}, or output of \code{\link{SimulateGraphical}},
#'   \code{\link{SimulateRegression}}.
#' @param graph2 second graph.
#' @param colours vector of edge colours. The first entry of the vector defines
#'   the colour of edges in \code{graph1} only, second entry is for edges in
#'   \code{graph2} only and third entry is for common edges.
#' @param lty vector of line types for edges. The order is defined as for
#'   argument \code{colours}.
#' @param show_labels logical indicating if the node labels should be displayed.
#' @param ... additional arguments to be passed to \code{\link{Graph}}.
#'
#' @seealso \code{\link{SelectionPerformanceGraph}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul1 <- SimulateGraphical(pk = 30)
#' set.seed(2)
#' simul2 <- SimulateGraphical(pk = 30)
#'
#' # Edge-wise comparison of the two graphs
#' mygraph <- GraphComparison(
#'   graph1 = simul1,
#'   graph2 = simul2
#' )
#' plot(mygraph, layout = layout_with_kk(mygraph))
#' }
#'
#' @export
GraphComparison <- function(graph1, graph2,
                            colours = c("tomato", "forestgreen", "navy"),
                            lty = c(2, 3, 1),
                            node_colour = NULL,
                            show_labels = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Extracting the adjacency matrices
  theta <- AdjacencyFromObject(graph1)
  theta_star <- AdjacencyFromObject(graph2)

  # Defining the vector of node colours
  if (is.null(node_colour)) {
    if (!"matrix" %in% class(graph1)) {
      if (class(graph1) == "variable_selection") {
        node_colour <- c(rep("skyblue", nrow(theta)), "red")
      }
      if (class(graph1) == "bi_selection") {
        if ("selected" %in% names(theta)) {
          selected <- theta$selected # PLS
        } else {
          selected <- theta$selectedX # PCA
        }
        node_colour <- c(
          rep("skyblue", nrow(selected)),
          rep("red", ncol(selected))
        )
      }
    }
    if (is.null(node_colour)) {
      node_colour <- "skyblue"
    }
  }

  # Checking input is a matrix
  if ((length(dim(theta)) != 2) | (length(dim(theta_star)) != 2)) {
    stop("Invalid input. Please provide adjacency matrices, igraph or sharp objects.")
  }

  # Extracting the estimated edges only
  if (ncol(theta) != ncol(theta_star)) {
    theta_star <- theta_star[rownames(theta), colnames(theta)]
  }

  # Storing similarities/differences between estimated and true sets
  Asum <- theta + 2 * theta_star

  # Refining inputs
  names(colours) <- names(lty) <- c("FP", "TN", "TP")

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sharp::Graph)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("adjacency", "node_colour")]

  # Making consensus graph
  # g <- Graph(adjacency = ifelse(Asum != 0, yes = 1, no = 0), node_colour = node_colour, ...)
  g <- do.call(sharp::Graph, args = c(
    list(
      adjacency = ifelse(Asum != 0, yes = 1, no = 0),
      node_colour = node_colour
    ),
    tmp_extra_args
  ))

  # Formatting vertices
  igraph::V(g)$size <- igraph::V(g)$size / 3 + 1
  if (!show_labels) {
    igraph::V(g)$label <- rep("", length(igraph::V(g)$label))
  }

  # Formatting edges
  myedgecolour <- colours[Asum[igraph::get.edgelist(g)]]
  myedgelty <- lty[Asum[igraph::get.edgelist(g)]]
  igraph::E(g)$color <- myedgecolour
  igraph::E(g)$width <- 1
  igraph::E(g)$lty <- myedgelty

  # Returning output graph
  return(g)
}


#' Graph representation of selection performance
#'
#' Generates an igraph object representing the True Positive, False Positive and
#' False Negative edges by comparing the set of selected edges to the set of
#' true edges. This function can only be used in simulation studies (i.e. when
#' the true model is known).
#'
#' @inheritParams GraphComparison
#' @param theta binary adjacency matrix or output of
#'   \code{\link{GraphicalModel}}, \code{\link{VariableSelection}}, or
#'   \code{\link{BiSelection}}.
#' @param theta_star true binary adjacency matrix or output of
#'   \code{\link{SimulateGraphical}} or \code{\link{SimulateRegression}}.
#' @param colours vector of edge colours. The first entry of the vector defines
#'   the colour of False Positive edges, second entry is for True Negatives and
#'   third entry is for True Positives.
#'
#' @family functions for model performance
#'
#' @seealso \code{\link{SimulateGraphical}}, \code{\link{SimulateRegression}},
#'   \code{\link{GraphicalModel}}, \code{\link{VariableSelection}},
#'   \code{\link{BiSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 30)
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data, K = 10)
#'
#' # Performance graph
#' perfgraph <- SelectionPerformanceGraph(
#'   theta = stab,
#'   theta_star = simul
#' )
#' plot(perfgraph)
#'
#' # Alternative formulation
#' perfgraph <- SelectionPerformanceGraph(
#'   theta = Adjacency(stab),
#'   theta_star = simul$theta
#' )
#' plot(perfgraph)
#'
#' # User-defined colours/shapes
#' perfgraph <- SelectionPerformanceGraph(
#'   theta = stab, theta_star = simul,
#'   colours = c("forestgreen", "orange", "black"),
#'   node_colour = "red", node_shape = "star"
#' )
#' plot(perfgraph)
#' perfgraph <- SelectionPerformanceGraph(
#'   theta = stab, theta_star = simul,
#'   colours = c("forestgreen", "orange", "black"), lty = c(4, 2, 3)
#' )
#' plot(perfgraph)
#'
#' # Using and re-formatting igraph object
#' require(igraph)
#' igraph::V(perfgraph)$size <- 10
#' plot(perfgraph, layout = igraph::layout_with_kk(perfgraph))
#'
#' # Regression model
#' set.seed(1)
#' simul <- SimulateRegression(pk = 30)
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' perf <- SelectionPerformance(theta = stab, theta_star = simul)
#' perf_graph <- SelectionPerformanceGraph(theta = stab, theta_star = simul)
#' plot(perf_graph)
#'
#' # Sparse PLS model
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#' stab <- BiSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(x) - 1),
#'   LambdaY = 1:(ncol(y) - 1),
#'   implementation = SparsePLS,
#'   n_cat = 2
#' )
#' perf <- SelectionPerformance(theta = stab, theta_star = simul)
#' perf_graph <- SelectionPerformanceGraph(theta = stab, theta_star = simul)
#' plot(perf_graph)
#' }
#'
#' @export
SelectionPerformanceGraph <- function(theta, theta_star,
                                      colours = c("tomato", "forestgreen", "navy"),
                                      lty = c(2, 3, 1),
                                      node_colour = NULL,
                                      show_labels = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  g <- GraphComparison(
    graph1 = theta,
    graph2 = theta_star,
    colours = colours,
    lty = lty,
    node_colour = node_colour,
    show_labels = show_labels,
    ...
  )
  #
  # # Re-formatting input theta
  # if (any(class(theta) %in% c("variable_selection", "graphical_model", "bi_selection"))) {
  #   if (class(theta) == "variable_selection") {
  #     theta <- cbind(SelectedVariables(theta))
  #     if (is.null(node_colour)) {
  #       node_colour <- c(rep("skyblue", nrow(theta)), "red")
  #     }
  #     if (ncol(theta) == 1) {
  #       colnames(theta) <- "outcome"
  #     }
  #     theta <- Square(theta)
  #   } else {
  #     if (class(theta) == "graphical_model") {
  #       theta <- Adjacency(theta)
  #     } else {
  #       if ("selected" %in% names(theta)) {
  #         selected <- theta$selected # PLS
  #       } else {
  #         selected <- theta$selectedX # PCA
  #       }
  #       # Relating predictors and outcomes (ignoring latent variables)
  #       if (is.null(node_colour)) {
  #         node_colour <- c(
  #           rep("skyblue", nrow(selected)),
  #           rep("red", ncol(selected))
  #         )
  #       }
  #       theta <- Square(selected)
  #     }
  #   }
  # }
  #
  # # Re-formatting input theta_star
  # if (any(class(theta_star) %in% c("simulation_regression", "simulation_graphical_model", "simulation_components"))) {
  #   if (class(theta_star) %in% c("simulation_regression", "simulation_components")) {
  #     theta_star <- cbind(theta_star$theta)
  #     if (ncol(theta_star) == 1) {
  #       colnames(theta_star) <- "outcome"
  #     }
  #     theta_star <- Square(theta_star)
  #   } else {
  #     theta_star <- theta_star$theta
  #   }
  # }
  #
  # # Defining node colour if not done before
  # if (is.null(node_colour)) {
  #   node_colour <- "skyblue"
  # }
  #
  # # Checking input is a matrix
  # if ((length(dim(theta)) != 2) | (length(dim(theta_star)) != 2)) {
  #   stop("Arguments 'theta' and 'theta_star' must be adjacency matrices.")
  # }
  #
  # # Extracting the estimated edges only
  # if (ncol(theta) != ncol(theta_star)) {
  #   theta_star <- theta_star[rownames(theta), colnames(theta)]
  # }
  #
  # # Storing similarities/differences between estimated and true sets
  # Asum <- theta + 2 * theta_star
  #
  # # Refining inputs
  # names(colours) <- names(lty) <- c("FP", "TN", "TP")
  #
  # # Extracting relevant extra arguments
  # tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = sharp::Graph)
  # tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("adjacency", "node_colour")]
  #
  # # Making consensus graph
  # # g <- Graph(adjacency = ifelse(Asum != 0, yes = 1, no = 0), node_colour = node_colour, ...)
  # g <- do.call(sharp::Graph, args = c(list(adjacency = ifelse(Asum != 0, yes = 1, no = 0), node_colour = node_colour), tmp_extra_args))
  #
  # # Formatting vertices
  # igraph::V(g)$size <- igraph::V(g)$size / 3 + 1
  # if (!show_labels) {
  #   igraph::V(g)$label <- rep("", length(igraph::V(g)$label))
  # }
  #
  # # Formatting edges
  # myedgecolour <- colours[Asum[igraph::get.edgelist(g)]]
  # myedgelty <- lty[Asum[igraph::get.edgelist(g)]]
  # igraph::E(g)$color <- myedgecolour
  # igraph::E(g)$width <- 1
  # igraph::E(g)$lty <- myedgelty

  # Returning output graph
  return(g)
}


#' Selection performance (internal)
#'
#' Computes different metrics of selection performance from a categorical
#' vector/matrix with 3 for True Positive, 2 for False Negative, 1 for False
#' Positive and 0 for True Negative.
#'
#' @inheritParams SelectionPerformance
#' @param Asum vector (in variable selection) or matrix (in graphical modelling)
#'   containing values of \code{0}, \code{1}, \code{2} or \code{3}.
#'
#' @return A matrix of selection metrics including:
#'
#'   \item{TP}{number of True Positives (TP)} \item{FN}{number of False
#'   Negatives (TN)} \item{FP}{number of False Positives (FP)} \item{TN}{number
#'   of True Negatives (TN)} \item{sensitivity}{sensitivity, i.e. TP/(TP+FN)}
#'   \item{specificity}{specificity, i.e. TN/(TN+FP)} \item{accuracy}{accuracy,
#'   i.e. (TP+TN)/(TP+TN+FP+FN)} \item{precision}{precision (p), i.e.
#'   TP/(TP+FP)} \item{recall}{recall (r), i.e. TP/(TP+FN)}
#'   \item{F1_score}{F1-score, i.e. 2*p*r/(p+r)}
#'
#'   If argument "cor" is provided, the number of False Positives among
#'   correlated (FP_c) and uncorrelated (FP_i) pairs, defined as having
#'   correlations (provided in "cor") above or below the threshold "thr", are
#'   also reported.
#'
#' @keywords internal
SelectionPerformanceSingle <- function(Asum, cor = NULL, thr = 0.5) {
  # Asum is an adjacency matrix with 3 for TP, 2 for FN, 1 for FP, and 0 for TN

  # Preparing objects
  if (is.matrix(Asum)) {
    p <- ncol(Asum)
    N <- p * (p - 1) / 2
    Asum <- Asum[upper.tri(Asum)]
  } else {
    N <- length(Asum)
  }

  # Computing the numbers of True/False Positives/Negatives
  TP <- sum(Asum == 3)
  FN <- sum(Asum == 2)
  FP <- sum(Asum == 1)
  TN <- sum(Asum == 0)

  # Separation between correlated and independent features based on a threshold in correlation
  if (!is.null(cor)) {
    if (is.matrix(cor)) {
      cor_vect <- cor[upper.tri(cor)]
    } else {
      cor_vect <- cor
    }
    FP_c <- sum((Asum == 1) & (abs(cor_vect) >= thr))
    FP_i <- sum((Asum == 1) & (abs(cor_vect) < thr))
  }

  # Computing performances in selection
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / N
  if (TP + FP > 0) {
    precision <- TP / (TP + FP)
  } else {
    precision <- 0
  }
  if ((TP + FN) > 0) {
    recall <- TP / (TP + FN)
  } else {
    recall <- 1
  }
  if ((precision > 0) | (recall > 0)) {
    F1_score <- 2 * precision * recall / (precision + recall)
  } else {
    F1_score <- 0
  }

  if (is.null(cor)) {
    return(data.frame(
      TP = TP, FN = FN, FP = FP, TN = TN,
      sensitivity = sensitivity, specificity = specificity,
      accuracy = accuracy, precision = precision, recall = recall, F1_score = F1_score
    ))
  } else {
    return(data.frame(
      TP = TP, FN = FN, FP = FP, TN = TN, FP_c = FP_c, FP_i = FP_i,
      sensitivity = sensitivity, specificity = specificity,
      accuracy = accuracy, precision = precision, recall = recall, F1_score = F1_score
    ))
  }
}

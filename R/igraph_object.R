#' Graph visualisation
#'
#' Produces an \code{\link[igraph:igraph-package]{igraph}} object from an
#' adjacency matrix.
#'
#' @param adjacency adjacency matrix or output of \code{\link{GraphicalModel}}.
#' @param node_label optional vector of node labels. This vector must contain as
#'   many entries as there are rows/columns in the adjacency matrix and must be
#'   in the same order (the order is used to assign labels to nodes).
#' @param node_colour optional vector of node colours. This vector must contain
#'   as many entries as there are rows/columns in the adjacency matrix and must
#'   be in the same order (the order is used to assign colours to nodes).
#'   Integers, named colours or RGB values can be used.
#' @param node_shape optional vector of node shapes. This vector must contain as
#'   many entries as there are rows/columns in the adjacency matrix and must be
#'   in the same order (the order is used to assign shapes to nodes). Possible
#'   values are \code{"circle"}, \code{"square"}, \code{"triangle"} or
#'   \code{"star"}.
#' @param edge_colour optional character string for edge colour. Integers, named
#'   colours or RGB values can be used.
#' @param label_colour optional character string for label colour. Integers,
#'   named colours or RGB values can be used.
#' @param mode character string indicating how the adjacency matrix should be
#'   interpreted. Possible values include \code{"undirected"} or
#'   \code{"directed"} (see \code{\link[igraph]{graph_from_adjacency_matrix}}).
#' @param weighted indicating if entries of the adjacency matrix should define
#'   edge width. If \code{weighted=FALSE}, an unweighted igraph object is
#'   created, all edges have the same width. If \code{weighted=TRUE}, edge width
#'   is defined by the corresponding value in the adjacency matrix. If
#'   \code{weighted=NULL}, nodes are linked by as many edges as indicated in the
#'   adjacency matrix (integer values are needed).
#' @param satellites logical indicating if unconnected nodes (satellites) should
#'   be included in the igraph object.
#'
#' @return An igraph object.
#'
#' @details All functionalities implemented in
#'   \code{\link[igraph:igraph-package]{igraph}} can be used on the output.
#'   These include cosmetic changes for the visualisation, but also various
#'   tools for network analysis (including topological properties and community
#'   detection).
#'
#'   The R package \code{\link[visNetwork:visDocumentation]{visNetwork}} offers
#'   interactive network visualisation tools. An
#'   \code{\link[igraph:igraph-package]{igraph}} object can easily be converted
#'   to a \code{\link[visNetwork:visDocumentation]{visNetwork}} object (see
#'   example below).
#'
#'   For Cytoscape users, the \code{\link[RCy3:RCy3]{RCy3}} package can be used
#'   to open the network in Cytoscape.
#'
#' @seealso \code{\link{Adjacency}}, \code{\link{GraphicalModel}},
#'   \href{https://igraph.org/r/}{igraph manual},
#'   \href{http://datastorm-open.github.io/visNetwork/}{visNetwork manual},
#'   \href{https://cytoscape.org}{Cytoscape}
#'
#' @examples
#' \donttest{
#' ## From adjacency matrix
#'
#' # Un-weighted
#' adjacency <- SimulateAdjacency(pk = 20, topology = "scale-free")
#' plot(Graph(adjacency))
#'
#' # Weighted
#' adjacency <- adjacency * runif(prod(dim(adjacency)))
#' adjacency <- adjacency + t(adjacency)
#' plot(Graph(adjacency, weighted = TRUE))
#'
#' # Node colours and shapes
#' plot(Graph(adjacency, weighted = TRUE, node_shape = "star", node_colour = "red"))
#'
#'
#' ## From stability selection outputs
#'
#' # Graphical model
#' set.seed(1)
#' simul <- SimulateGraphical(pk = 20)
#' stab <- GraphicalModel(xdata = simul$data)
#' plot(Graph(stab))
#'
#' # Sparse PLS
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 50, pk = c(5, 5, 5), family = "gaussian")
#'   x <- simul$xdata
#'   y <- simul$ydata
#'   stab <- BiSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     family = "gaussian", ncomp = 3,
#'     LambdaX = seq_len(ncol(x) - 1),
#'     implementation = SparsePLS
#'   )
#'   plot(Graph(stab))
#' }
#'
#'
#' ## Tools from other packages
#'
#' # Applying some igraph functionalities
#' adjacency <- SimulateAdjacency(pk = 20, topology = "scale-free")
#' mygraph <- Graph(adjacency)
#' igraph::degree(mygraph)
#' igraph::betweenness(mygraph)
#' igraph::shortest_paths(mygraph, from = 1, to = 2)
#' igraph::walktrap.community(mygraph)
#'
#' # Interactive view using visNetwork
#' if (requireNamespace("visNetwork", quietly = TRUE)) {
#'   vgraph <- mygraph
#'   igraph::V(vgraph)$shape <- rep("dot", length(igraph::V(vgraph)))
#'   v <- visNetwork::visIgraph(vgraph)
#'   mylayout <- as.matrix(v$x$nodes[, c("x", "y")])
#'   mylayout[, 2] <- -mylayout[, 2]
#'   plot(mygraph, layout = mylayout)
#' }
#'
#' # Opening in Cytoscape using RCy3
#' if (requireNamespace("RCy3", quietly = TRUE)) {
#'   # Make sure that Cytoscape is open before running the following line
#'   # RCy3::createNetworkFromIgraph(mygraph)
#' }
#' }
#'
#' @export
Graph <- function(adjacency, node_label = NULL, node_colour = NULL, node_shape = NULL,
                  edge_colour = "grey60", label_colour = "grey20",
                  mode = "undirected", weighted = FALSE, satellites = FALSE) {
  # Checking input values (weighted)
  if (!is.null(weighted)) {
    if (!weighted %in% c(TRUE, FALSE)) {
      stop("Invalid input for argument 'weighted'. Possible values are: NULL, TRUE or FALSE.")
    }
  }

  # Extracting the adjacency matrix from the output of GraphicalModel()
  if (!is.matrix(adjacency)) {
    stability <- adjacency
    adjacency <- Adjacency(stability = adjacency)
  }

  # Setting row and column names if none
  if (is.null(rownames(adjacency)) & (is.null(colnames(adjacency)))) {
    rownames(adjacency) <- colnames(adjacency) <- paste0("var", seq_len(ncol(adjacency)))
  } else {
    if (is.null(rownames(adjacency))) {
      rownames(adjacency) <- colnames(adjacency)
    }
    if (is.null(colnames(adjacency))) {
      colnames(adjacency) <- rownames(adjacency)
    }
  }

  # Checking input values (node label)
  if (!is.null(node_label)) {
    if (length(node_label) != ncol(adjacency)) {
      stop(paste0(
        "Invalid input for argument 'node_label'. It must be a vector of length ",
        ncol(adjacency), " (the same as the number of nodes in the adjacency matrix)."
      ))
    }
  }

  # Checking input values (node colour)
  if (!is.null(node_colour)) {
    if (length(node_colour) == 1) {
      node_colour <- rep(node_colour, ncol(adjacency))
    } else {
      if (length(node_colour) != ncol(adjacency)) {
        stop(paste0(
          "Invalid input for argument 'node_colour'. It must be a vector of length ",
          ncol(adjacency), " (the same as the number of nodes in the adjacency matrix)."
        ))
      }
    }
  }

  # Checking input values (node shape)
  if (!is.null(node_shape)) {
    if (length(node_shape) == 1) {
      node_shape <- rep(node_shape, ncol(adjacency))
    } else {
      if (length(node_shape) != ncol(adjacency)) {
        stop(paste0(
          "Invalid input for argument 'node_shape'. It must be a vector of length ",
          ncol(adjacency), " (the same as the number of nodes in the adjacency matrix)."
        ))
      }
    }

    if (!any(node_shape %in% c("circle", "square", "triangle", "star"))) {
      stop(paste0("Invalid input for argument 'node_shape'. Possible values for the entries of the vector are: circle, square, triangle or star."))
    }
  }

  # Adding shapes if required
  if (!is.null(node_shape)) {
    if (any(node_shape == "star")) {
      igraph::add_shape("star",
        clip = igraph::shape_noclip,
        plot = mystar, parameters = list(vertex.norays = 5)
      )
    }

    if (any(node_shape == "triangle")) {
      igraph::add_shape("triangle",
        clip = igraph::shape_noclip,
        plot = mytriangle
      )
    }
  }

  # Default node colours
  if (is.null(node_colour)) {
    node_colour <- rep("skyblue", ncol(adjacency))
    if (exists("stability")) {
      if (inherits(stability, "bi_selection")) {
        node_colour <- ifelse(grepl("comp", colnames(adjacency)),
          yes = "orange", no = "skyblue"
        )
        node_colour[colnames(adjacency) %in% rownames(stability$selectedY)] <- "red"
      }
    }
  }

  # Default node shapes
  if (is.null(node_shape)) {
    node_shape <- rep("circle", ncol(adjacency))
  }

  # Default node labels
  if (is.null(node_label)) {
    node_label <- colnames(adjacency)
  }

  # Formatting node characteristics
  names(node_colour) <- colnames(adjacency)
  names(node_label) <- colnames(adjacency)
  names(node_shape) <- colnames(adjacency)

  # Formatting adjacency matrix
  if (!is.null(weighted)) {
    if (!weighted) {
      adjacency <- ifelse(adjacency != 0, yes = 1, no = 0)
      weighted <- NULL
    }
  } else {
    adjacency <- round(adjacency)
  }

  # Estimating igraph object
  mygraph <- igraph::graph_from_adjacency_matrix(adjacency, mode = mode, weighted = weighted)
  mydegrees <- igraph::degree(mygraph)

  # Changing arrow size for directed graphs
  if (mode == "directed") {
    igraph::E(mygraph)$arrow.size <- 1
  }

  # Including/excluding satellites (nodes with no edges)
  if (!satellites) {
    mygraph <- igraph::delete.vertices(mygraph, v = names(mydegrees)[mydegrees == 0])
  }

  # Formatting vertices
  mydegrees <- igraph::degree(mygraph)
  igraph::V(mygraph)$size <- 10
  igraph::V(mygraph)$label <- node_label[igraph::V(mygraph)$name]
  igraph::V(mygraph)$color <- node_colour[igraph::V(mygraph)$name]
  igraph::V(mygraph)$shape <- node_shape[igraph::V(mygraph)$name]
  igraph::V(mygraph)$frame.color <- igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family <- "sans"
  igraph::V(mygraph)$label.cex <- 1
  igraph::V(mygraph)$label.color <- label_colour

  # Formatting edges
  igraph::E(mygraph)$color <- edge_colour
  if (is.null(weighted)) {
    igraph::E(mygraph)$width <- 0.5
  } else {
    igraph::E(mygraph)$width <- igraph::E(mygraph)$weight
  }

  return(mygraph)
}


#' Star-shaped nodes
#'
#' Produces star-shaped nodes in an igraph object.
#'
#' @param coords a matrix of coordinates
#' (see \code{\link[igraph]{add_shape}}).
#' @param v a vector of node IDs
#' (see \code{\link[igraph]{add_shape}}).
#' @param params node graphical parameters
#' (see \code{\link[igraph]{add_shape}}).
#'
#' @seealso \code{\link[igraph]{add_shape}}
#'
#' @keywords internal
mystar <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if ((length(vertex.color) != 1) & !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if ((length(vertex.size) != 1) & !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if ((length(norays) != 1) & !is.null(v)) {
    norays <- norays[v]
  }

  mapply(coords[, 1], coords[, 2], vertex.color, vertex.size, norays,
    FUN = function(x, y, bg, size, nor) {
      graphics::symbols(
        x = x, y = y, fg = bg, bg = bg,
        stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
        add = TRUE, inches = FALSE
      )
    }
  )
}


#' Triangular nodes
#'
#' Produces triangular nodes in an igraph object.
#'
#' @param coords a matrix of coordinates
#' (see \code{\link[igraph]{add_shape}}).
#' @param v a vector of node IDs
#' (see \code{\link[igraph]{add_shape}}).
#' @param params node graphical parameters
#' (see \code{\link[igraph]{add_shape}}).
#'
#' @seealso \code{\link[igraph]{add_shape}}
#'
#' @keywords internal
mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if ((length(vertex.color) != 1) & !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if ((length(vertex.size) != 1) & !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- 3

  mapply(coords[, 1], coords[, 2], vertex.color, vertex.size, norays,
    FUN = function(x, y, bg, size, nor) {
      graphics::symbols(
        x = x, y = y, fg = bg, bg = bg,
        stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
        add = TRUE, inches = FALSE
      )
    }
  )
}

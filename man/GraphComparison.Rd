% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_performance.R
\name{GraphComparison}
\alias{GraphComparison}
\title{Edge-wise comparison of two graphs}
\usage{
GraphComparison(
  graph1,
  graph2,
  col = c("tomato", "forestgreen", "navy"),
  lty = c(2, 3, 1),
  node_colour = NULL,
  show_labels = TRUE,
  ...
)
}
\arguments{
\item{graph1}{first graph. Possible inputs are: adjacency matrix, or
\code{\link[igraph:igraph-package]{igraph}} object, or output of
\code{\link{GraphicalModel}}, \code{\link{VariableSelection}},
\code{\link{BiSelection}}, or output of \code{\link[fake]{SimulateGraphical}},
\code{\link[fake]{SimulateRegression}}.}

\item{graph2}{second graph.}

\item{col}{vector of edge colours. The first entry of the vector defines
the colour of edges in \code{graph1} only, second entry is for edges in
\code{graph2} only and third entry is for common edges.}

\item{lty}{vector of line types for edges. The order is defined as for
argument \code{col}.}

\item{node_colour}{optional vector of node colours. This vector must contain
as many entries as there are rows/columns in the adjacency matrix and must
be in the same order (the order is used to assign colours to nodes).
Integers, named colours or RGB values can be used.}

\item{show_labels}{logical indicating if the node labels should be displayed.}

\item{...}{additional arguments to be passed to \code{\link{Graph}}.}
}
\value{
An igraph object.
}
\description{
Generates an \code{\link[igraph:igraph-package]{igraph}} object representing
the common and graph-specific edges.
}
\examples{
# Data simulation
set.seed(1)
simul1 <- SimulateGraphical(pk = 30)
set.seed(2)
simul2 <- SimulateGraphical(pk = 30)

# Edge-wise comparison of the two graphs
mygraph <- GraphComparison(
  graph1 = simul1,
  graph2 = simul2
)
plot(mygraph, layout = igraph::layout_with_kk(mygraph))
}
\seealso{
\code{\link{SelectionPerformanceGraph}}
}

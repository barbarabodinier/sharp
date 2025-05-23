% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_performance.R
\name{SelectionPerformanceGraph}
\alias{SelectionPerformanceGraph}
\title{Graph representation of selection performance}
\usage{
SelectionPerformanceGraph(
  theta,
  theta_star,
  col = c("tomato", "forestgreen", "navy"),
  lty = c(2, 3, 1),
  node_colour = NULL,
  show_labels = TRUE,
  ...
)
}
\arguments{
\item{theta}{binary adjacency matrix or output of
\code{\link{GraphicalModel}}, \code{\link{VariableSelection}}, or
\code{\link{BiSelection}}.}

\item{theta_star}{true binary adjacency matrix or output of
\code{\link[fake]{SimulateGraphical}} or \code{\link[fake]{SimulateRegression}}.}

\item{col}{vector of edge colours. The first entry of the vector defines
the colour of False Positive edges, second entry is for True Negatives and
third entry is for True Positives.}

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
Generates an igraph object representing the True Positive, False Positive and
False Negative edges by comparing the set of selected edges to the set of
true edges. This function can only be used in simulation studies (i.e. when
the true model is known).
}
\examples{
\donttest{
# Data simulation
set.seed(1)
simul <- SimulateGraphical(pk = 30)

# Stability selection
stab <- GraphicalModel(xdata = simul$data, K = 10)

# Performance graph
perfgraph <- SelectionPerformanceGraph(
  theta = stab,
  theta_star = simul
)
plot(perfgraph)
}

}
\seealso{
\code{\link[fake]{SimulateGraphical}}, \code{\link[fake]{SimulateRegression}},
  \code{\link{GraphicalModel}}, \code{\link{VariableSelection}},
  \code{\link{BiSelection}}

Other functions for model performance: 
\code{\link{ClusteringPerformance}()},
\code{\link{SelectionPerformance}()}
}
\concept{functions for model performance}

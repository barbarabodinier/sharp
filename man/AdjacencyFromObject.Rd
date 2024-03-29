% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{AdjacencyFromObject}
\alias{AdjacencyFromObject}
\title{Adjacency matrix from object}
\usage{
AdjacencyFromObject(object)
}
\arguments{
\item{object}{input object.}
}
\value{
A vector without missing values or NULL.
}
\description{
Returns the adjacency matrix from an
\code{\link[igraph:igraph-package]{igraph}} object or from the output of
simulation or inference functions from the present package.
}
\keyword{internal}

#' Simulation of data with underlying graph structure
#'
#' Simulates (i) a graph, and (ii) multivariate Normal data for which the graph
#' structure is encoded in the nonzero entries of the true partial correlation
#' matrix. This procedure ensures that the conditional independence structure
#' between the variables is encoded in the simulated graph. The outputs of this
#' function can be used to evaluate the ability of a graphical model to identify
#' edges of a conditional independence graph.
#'
#' @param n number of observations in the simulated data.
#' @param pk vector of the number of variables per group in the simulated data.
#'   The number of nodes in the simulated graph is \code{sum(pk)}. With multiple
#'   groups, the simulated (partial) correlation matrix has a block structure,
#'   where blocks arise from the integration of the \code{length(pk)} groups.
#'   This argument is only used if \code{sum(pk)} is equal to the number of
#'   rows/columns in \code{theta} is not provided.
#' @param theta optional binary and symmetric adjacency matrix encoding the
#'   conditional independence structure.
#' @param implementation function for simulation of the graph. By default,
#'   algorithms implemented in \code{\link[huge]{huge.generator}} are used.
#'   Alternatively, a user-defined function can be used. It must take \code{pk},
#'   \code{topology} and \code{nu} as arguments and return a
#'   \code{(sum(pk)*(sum(pk)))} binary and symmetric matrix for which diagonal
#'   entries are all equal to zero. This function is only applied if
#'   \code{theta} is not provided.
#' @param topology topology of the simulated graph. If using
#'   \code{implementation=HugeAdjacency}, possible values are listed for the
#'   argument \code{graph} of \code{\link[huge]{huge.generator}}. These are:
#'   "random", "hub", "cluster", "band" and "scale-free".
#' @param nu_within expected density (number of edges over the number of node
#'   pairs) of within-group blocks in the graph. If \code{length(pk)=1}, this is
#'   the expected density of the graph. If \code{implementation=HugeAdjacency},
#'   this argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}).
#' @param nu_between expected density (number of edges over the number of node
#'   pairs) of between-group blocks in the graph. Similar to \code{nu_within}.
#'   By default, the same density is used for within and between blocks
#'   (\code{nu_within}=\code{nu_between}). Only used if \code{length(pk)>1}.
#' @param output_matrices logical indicating if the true precision and (partial)
#'   correlation matrices should be included in the output.
#' @param v_within vector defining the (range of) nonzero entries in the
#'   diagonal blocks of the precision matrix. These values must be between -1
#'   and 1 if \code{pd_strategy="min_eigenvalue"}. If \code{continuous=FALSE},
#'   \code{v_within} is the set of possible precision values. If
#'   \code{continuous=TRUE}, \code{v_within} is the range of possible precision
#'   values.
#' @param v_between vector defining the (range of) nonzero entries in the
#'   off-diagonal blocks of the precision matrix. This argument is the same as
#'   \code{v_within} but for off-diagonal blocks. It is only used if
#'   \code{length(pk)>1}.
#' @param v_sign vector of possible signs for precision matrix entries. Possible
#'   inputs are: \code{-1} for positive partial correlations, \code{1} for
#'   negative partial correlations, or \code{c(-1, 1)} for both positive and
#'   negative partial correlations.
#' @param continuous logical indicating whether to sample precision values from
#'   a uniform distribution between the minimum and maximum values in
#'   \code{v_within} (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=TRUE}) or from proposed values in \code{v_within}
#'   (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=FALSE}).
#' @param pd_strategy method to ensure that the generated precision matrix is
#'   positive definite (and hence can be a covariance matrix). With
#'   \code{pd_strategy="diagonally_dominant"}, the precision matrix is made
#'   diagonally dominant by setting the diagonal entries to the sum of absolute
#'   values on the corresponding row and a constant u. With
#'   \code{pd_strategy="min_eigenvalue"}, diagonal entries are set to the sum of
#'   the absolute value of the smallest eigenvalue of the precision matrix with
#'   zeros on the diagonal and a constant u.
#' @param ev expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis. This is the largest
#'   eigenvalue of the correlation (if \code{scale=TRUE}) or covariance (if
#'   \code{scale=FALSE}) matrix divided by the sum of eigenvalues. If
#'   \code{ev=NULL} (the default), the constant u is chosen by maximising the
#'   contrast of the correlation matrix.
#' @param scale logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix. If \code{scale=TRUE}, the
#'   correlation matrix is used as parameter of the multivariate normal
#'   distribution.
#' @param u_list vector with two numeric values defining the range of values to
#'   explore for constant u.
#' @param tol accuracy for the search of parameter u as defined in
#'   \code{\link[stats]{optimise}}.
#' @param ... additional arguments passed to the graph simulation function
#'   provided in \code{implementation}.
#'
#' @seealso \code{\link{SimulatePrecision}}, \code{\link{MakePositiveDefinite}},
#'   \code{\link{Contrast}}, \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables.} \item{theta}{adjacency matrix of the simulated
#'   graph} \item{omega}{simulated (true) precision matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{phi}{simulated (true) partial
#'   correlation matrix. Only returned if \code{output_matrices=TRUE}.}
#'   \item{sigma}{ simulated (true) covariance matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u}{value of
#'   the constant u used for the simulation of \code{omega}. Only returned if
#'   \code{output_matrices=TRUE}.}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of random graph with 50 nodes
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 50, topology = "random", nu_within = 0.05)
#' dim(simul$data) # dimension of simulated dataset
#' sum(simul$theta) / 2 # number of edges
#' plot(Graph(simul$theta))
#'
#' # Simulation of scale-free graph with 20 nodes
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, topology = "scale-free")
#' plot(Graph(simul$theta))
#'
#' # Extracting true precision/correlation matrices
#' set.seed(1)
#' simul <- SimulateGraphical(
#'   n = 100, pk = 20,
#'   topology = "scale-free", output_matrices = TRUE
#' )
#' str(simul)
#'
#' # Simulation of multi-block data
#' set.seed(1)
#' pk <- c(20, 30)
#' simul <- SimulateGraphical(
#'   n = 100, pk = pk,
#'   nu_within = 0.05, nu_between = 0.05,
#'   pd_strategy = "min_eigenvalue"
#' )
#' mycor <- cor(simul$data)
#' Heatmap(mycor,
#'   colours = c("darkblue", "white", "firebrick3"),
#'   legend_range = c(-1, 1), legend_length = 50,
#'   legend = FALSE, axes = FALSE
#' )
#' for (i in 1:2) {
#'   axis(side = i, at = c(0.5, pk[1] - 0.5), labels = NA)
#'   axis(
#'     side = i, at = mean(c(0.5, pk[1] - 0.5)),
#'     labels = ifelse(i == 1, yes = "Group 1", no = "Group 2"),
#'     tick = FALSE, cex.axis = 1.5
#'   )
#'   axis(side = i, at = c(pk[1] + 0.5, sum(pk) - 0.5), labels = NA)
#'   axis(
#'     side = i, at = mean(c(pk[1] + 0.5, sum(pk) - 0.5)),
#'     labels = ifelse(i == 1, yes = "Group 2", no = "Group 1"),
#'     tick = FALSE, cex.axis = 1.5
#'   )
#' }
#'
#' # User-defined function for graph simulation
#' CentralNode <- function(pk, hub = 1) {
#'   theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
#'   theta[hub, ] <- 1
#'   theta[, hub] <- 1
#'   diag(theta) <- 0
#'   return(theta)
#' }
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = CentralNode)
#' plot(Graph(simul$theta)) # star
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = CentralNode, hub = 2)
#' plot(Graph(simul$theta)) # variable 2 is the central node
#'
#' # User-defined adjacency matrix
#' mytheta <- matrix(c(
#'   0, 1, 1, 0,
#'   1, 0, 0, 0,
#'   1, 0, 0, 1,
#'   0, 0, 1, 0
#' ), ncol = 4, byrow = TRUE)
#' simul <- SimulateGraphical(n = 100, theta = mytheta)
#' plot(Graph(simul$theta))
#'
#' # User-defined adjacency and block structure
#' simul <- SimulateGraphical(n = 100, theta = mytheta, pk = c(2, 2))
#' mycor <- cor(simul$data)
#' Heatmap(mycor,
#'   colours = c("darkblue", "white", "firebrick3"),
#'   legend_range = c(-1, 1), legend_length = 50, legend = FALSE
#' )
#' }
#' @export
SimulateGraphical <- function(n = 100, pk = 10, theta = NULL,
                              implementation = HugeAdjacency, topology = "random",
                              nu_within = 0.1, nu_between = NULL,
                              v_within = c(0.5, 1), v_between = c(0, 0.1),
                              v_sign = c(-1, 1), continuous = TRUE,
                              pd_strategy = "diagonally_dominant", ev = NULL, scale = TRUE,
                              u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                              output_matrices = FALSE, ...) {
  # Defining number of nodes
  p <- sum(pk)
  if (!is.null(theta)) {
    if (ncol(theta) != p) {
      p <- pk <- ncol(theta)
    }
  }

  # Defining the between-block density
  if (is.null(nu_between)) {
    nu_between <- nu_within
  }

  # Building adjacency matrix
  if (is.null(theta)) {
    theta <- SimulateAdjacency(
      pk = pk,
      implementation = implementation, topology = topology,
      nu_within = nu_within, nu_between = nu_between, ...
    )
  }

  # Simulation of a precision matrix
  out <- SimulatePrecision(
    pk = pk, theta = theta,
    v_within = v_within, v_between = v_between,
    v_sign = v_sign, continuous = continuous,
    pd_strategy = pd_strategy, ev = ev, scale = scale,
    u_list = u_list, tol = tol
  )
  omega <- out$omega

  # Computing the covariance matrix
  if (scale) {
    sigma <- stats::cov2cor(solve(omega))
  } else {
    sigma <- solve(omega)
  }

  # Computing the partial correlation matrix
  if (output_matrices) {
    phi <- -stats::cov2cor(omega) + 2 * diag(ncol(omega))
  }

  # Simulating data from multivariate normal distribution
  x <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(x) <- paste0("var", 1:ncol(x))
  rownames(x) <- paste0("obs", 1:nrow(x))

  if (output_matrices) {
    return(list(
      data = x, theta = theta,
      omega = omega, phi = phi, sigma = sigma,
      u = out$u
    ))
  } else {
    return(list(data = x, theta = theta))
  }
}


#' Simulation of data with underlying clusters
#'
#' Simulates mixture multivariate Normal data with clusters of observations
#' (rows) sharing similar profiles along (a subset of) variables (columns). The
#' conditional independence structure between the variables can be simulated or
#' provided in argument \code{adjacency}. The same covariance is used across all
#' clusters. Independent variables are simulated by default
#' (\code{nu_within=0}).
#'
#' @inheritParams SimulateGraphical
#' @param n vector of the number of observations per cluster in the simulated
#'   data. The number of observations in the simulated data is \code{sum(n)}.
#' @param pk vector of the number of variables in the simulated data.
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional independence structure between variables.
#' @param theta_xc optional binary vector encoding which variables (columns)
#'   contribute to the clustering structure between observations (rows).
#' @param nu_xc expected proportion of variables contributing to the clustering
#'   over the total number of variables. This argument is only used if
#'   \code{theta_xc} is not provided.
#' @param ev_xc vector of marginal expected proportion of explained for each
#'   variable contributing to the clustering. This parameter is only used for
#'   variables with a nonzero entry in \code{theta_xc}.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis applied on the
#'   predictors. This is the largest eigenvalue of the correlation (if
#'   \code{scale=TRUE}) or covariance (if \code{scale=FALSE}) matrix divided by
#'   the sum of eigenvalues. If \code{ev=NULL} (the default), the constant u is
#'   chosen by maximising the contrast of the correlation matrix.
#'
#' @seealso \code{\link{MakePositiveDefinite}}, \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{sum(n)}
#'   observation and \code{sum(pk)} variables} \item{theta}{simulated (true)
#'   cluster membership.} \item{theta}{adjacency matrix of the graph encoding
#'   the conditional independence structure between variables.}
#'   \item{theta_xc}{binary vector encoding variables contributing to the
#'   clustering structure.} \item{ev}{vector of marginal expected proportions of
#'   explained variance for each variable.}
#'
#' @examples
#' \dontrun{
#' ## Example with 3 clusters
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(10, 30, 15)
#' )
#'
#' # Visualisation of Euclidian distances
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = as.matrix(dist(simul$data)),
#'   colours = c("navy", "white", "red")
#' )
#'
#'
#' ## Example with 2 variables contributing to clustering
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8))
#' )
#'
#' # Visualisation of the data
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = simul$data,
#'   colours = c("navy", "white", "red")
#' )
#' simul$ev # marginal proportions of explained variance
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with more distinct clusters
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = 10,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.9, 0.8, rep(0, 8))
#' )
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#'
#' ## Example with correlated contributors
#'
#' # Data simulation
#' pk <- 10
#' adjacency <- matrix(0, pk, pk)
#' adjacency[1, 2] <- adjacency[2, 1] <- 1
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(200, 100, 150), pk = pk,
#'   theta_xc = c(1, 1, rep(0, 8)),
#'   ev_xc = c(0.9, 0.8, rep(0, 8)),
#'   adjacency = adjacency,
#'   pd_strategy = "min_eigenvalue",
#'   v_within = 0.6, v_sign = -1
#' )
#'
#' # Visualisation along contributing variables
#' plot(simul$data[, 1:2], col = simul$theta)
#'
#' # Checking marginal proportions of explained variance
#' mymodel <- lm(simul$data[, 1] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#' mymodel <- lm(simul$data[, 2] ~ as.factor(simul$theta))
#' summary(mymodel)$r.squared
#' }
#' @export
SimulateClustering <- function(n = c(10, 10), pk = 10, adjacency = NULL,
                               theta_xc = NULL, nu_xc = 0.1, ev_xc = NULL,
                               implementation = HugeAdjacency, topology = "random",
                               nu_within = 0, nu_between = NULL,
                               v_within = c(0.5, 1), v_between = c(0, 0.1),
                               v_sign = c(-1, 1), continuous = TRUE,
                               pd_strategy = "diagonally_dominant", ev_xx = NULL, scale = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                               output_matrices = FALSE) {
  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = sum(n), pk = pk, theta = adjacency,
    implementation = implementation,
    topology = topology,
    nu_within = nu_within,
    nu_between = nu_between,
    output_matrices = output_matrices,
    v_within = v_within,
    v_between = v_between,
    continuous = continuous,
    pd_strategy = pd_strategy, ev = ev_xx, scale = scale,
    u_list = u_list, tol = tol
  )

  # Defining number of clusters
  nc <- length(n)

  # Defining variables contributing to the clustering
  if (is.null(theta_xc)) {
    theta_xc <- SamplePredictors(pk = sum(pk), q = 1, nu = nu_xc, orthogonal = TRUE)[, 1]
  }

  # Simulating marginal proportions of explained variance
  if (is.null(ev_xc)) {
    ev_xc <- stats::runif(n = sum(pk))
  } else {
    if (length(ev_xc) == 1) {
      ev_xc <- rep(ev_xc, sum(pk))
    }
  }

  # Re-naming the outputs
  out$adjacency <- out$theta

  # Definition of membership
  theta <- NULL
  for (i in 1:length(n)) {
    theta <- c(theta, rep(i, each = n[i]))
  }
  names(theta) <- rownames(out$data)
  out$theta <- theta

  # Building binary cluster membership for each feature
  V <- stats::model.matrix(~ as.factor(theta) - 1)

  # Simulating the cluster-specific means
  mu_mat <- matrix(NA, nrow = sum(n), ncol = sum(pk))
  for (k in 1:ncol(mu_mat)) {
    # Defining variance to reach expected proportion of e.v.
    var_mu <- ev_xc[k] * 1 / (1 - ev_xc[k])

    # Sampling initial values for cluster-specific means
    mu <- stats::rnorm(n = nc, mean = 0, sd = 1)
    for (i in 1:nrow(mu_mat)) {
      mu_mat[i, k] <- mu[theta[i]]
    }

    # Scaling to ensure mean of zero and defined variance
    mu_mat[, k] <- scale(mu_mat[, k])
    mu_mat[, k] <- mu_mat[, k] * sqrt(var_mu)
  }

  # Using cluster-specific mean for contributing variables
  # TODO: check impact of using non-scaled data here?
  for (k in 1:ncol(mu_mat)) {
    if (theta_xc[k] == 1) {
      out$data[, k] <- out$data[, k] + mu_mat[, k]
    }
  }

  # Definition of contributing variables
  names(theta_xc) <- colnames(out$data)
  out$theta_xc <- theta_xc
  out$ev <- ev_xc * theta_xc

  return(out)
}


#' Simulation of sparse orthogonal components
#'
#' Simulates variables following a multivariate Normal distribution that could
#' be obtained from a sparse linear combination of orthogonal latent variables.
#' This generates blocks of mutually independent variables, where all variables
#' from a block can be obtained from a linear combination of the same latent
#' variables. The latent variables would correspond to Principal Components from
#' a sparse Principal Component Analysis. The loadings coefficients, their
#' support, and the proportions of explained variance by each of the latent
#' variables are returned. This function can be used to evaluate the performance
#' of sparse Principal Component Analysis algorithms.
#'
#' @inheritParams SimulateGraphical
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional graph structure between observations. The clusters encoded in
#'   this argument must be in line with those indicated in \code{pk}. Edges in
#'   off-diagonal blocks are not allowed to ensure that the simulated orthogonal
#'   components are sparse. Corresponding entries in the precision matrix will
#'   be set to zero.
#'
#' @details The data is simulated from a centered multivariate Normal
#'   distribution with a block-diagonal covariance matrix. Independence between
#'   variables from the different blocks ensures that sparse orthogonal
#'   components can be generated. The block-diagonal (partial) correlation
#'   matrix is obtained using a graph structure encoding the conditional
#'   independence between variables. The orthogonal latent variables are
#'   obtained from eigendecomposition of the true correlation matrix. The sparse
#'   eigenvectors contain the weights of the linear combination of variables to
#'   construct the latent variable (loadings coefficients). The proportion of
#'   explained variance by each of the latent variable is computed from
#'   eigenvalues. As latent variables are defined from the true correlation
#'   matrix, the number of sparse orthogonal components is not limited by the
#'   number of observations and is equal to \code{sum(pk)}.
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables.} \item{loadings}{loadings coefficients of the
#'   orthogonal latent variables (principal components).} \item{theta}{support
#'   of the loadings coefficients.} \item{ev}{proportion of explained variance
#'   by each of the orthogonal latent variables.} \item{adjacency}{adjacency
#'   matrix of the simulated graph.} \item{omega}{simulated (true) precision
#'   matrix. Only returned if \code{output_matrices=TRUE}.} \item{phi}{simulated
#'   (true) partial correlation matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{C}{ simulated (true) correlation
#'   matrix. Only returned if \code{output_matrices=TRUE}.}
#'
#' @seealso \code{\link{MakePositiveDefinite}}, \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @examples
#' \dontrun{
#' # Simulation of 3 components with high e.v.
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4), ev = 0.4)
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(simul$data),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#' print(simul$ev)
#'
#' # Simulation of 3 components with moderate e.v.
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4), ev = 0.25)
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(simul$data),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#' print(simul$ev)
#'
#' # Simulation of multiple components with low e.v.
#' pk <- sample(3:10, size = 5, replace = TRUE)
#' simul <- SimulateComponents(
#'   pk = pk,
#'   nu_within = 0.3, v_within = c(0.8, 0.5), v_sign = -1, ev = 0.1
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(simul$data),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' plot(cumsum(simul$ev), ylim = c(0, 1), las = 1)
#' }
#' @export
SimulateComponents <- function(n = 100, pk = c(10, 10), adjacency = NULL,
                               nu_within = 1,
                               v_within = c(0.5, 1), v_sign = -1, continuous = TRUE,
                               pd_strategy = "min_eigenvalue", ev = 0.1, scale = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25,
                               output_matrices = FALSE) {
  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = n, pk = pk, theta = adjacency,
    implementation = HugeAdjacency,
    topology = "random",
    nu_within = nu_within, # fully connected components by default
    nu_between = 0, # need unconnected blocks
    v_within = v_within,
    v_between = 0,
    v_sign = v_sign,
    continuous = continuous,
    pd_strategy = pd_strategy, ev = ev, scale = scale,
    u_list = u_list, tol = tol,
    output_matrices = TRUE
  )

  # Eigendecomposition of the covariance
  eig <- eigen(out$sigma)

  # Definition of membership
  membership <- NULL
  for (i in 1:length(pk)) {
    membership <- c(membership, rep(i, each = pk[i]))
  }
  names(membership) <- colnames(out$data)
  out$membership <- membership

  # Re-naming the outputs
  out$adjacency <- out$theta

  # Definition of sparse principal components
  out$loadings <- round(eig$vectors, digits = 10)
  out$theta <- ifelse(out$loadings != 0, yes = 1, no = 0)
  rownames(out$theta) <- rownames(out$loadings) <- colnames(out$adjacency)
  colnames(out$theta) <- colnames(out$loadings) <- paste0("PC", 1:ncol(out$theta))

  # Definition of proportion of explained variance
  ev <- eig$values / sum(eig$values)
  names(ev) <- colnames(out$theta)
  out$ev <- ev

  # Re-arranging the output
  out <- out[c("data", "loadings", "theta", "ev", "membership", "omega", "phi", "C", "u")]
  if (!output_matrices) {
    out <- out[c("data", "loadings", "theta", "ev", "membership")]
  }

  return(out)
}


#' Simulation of predictors and associated outcome
#'
#' Simulates (i) a matrix \code{xdata} of \code{n} observations for
#' \code{sum(pk)} normally distributed predictor variables, (ii) a matrix
#' \code{zdata} of \code{length(pk)} orthogonal latent variables, and (iii) a
#' matrix \code{ydata} of \code{length(pk)} outcome variables. The conditional
#' independence structure between the predictors and latent variables is encoded
#' in a precision matrix, where the diagonal entries corresponding to latent
#' variables are tuned to reach a user-defined expected proportion of explained
#' variance. To ensure that latent variables are orthogonal (these can be
#' interpreted as the Principal Components of a Partial Least Squares model),
#' the predictors contributing to their definition are taken from independent
#' blocks of variables. The outcome variables are then obtained from a linear
#' combination of the latent variables. The outputs of this function can be used
#' to evaluate the ability of variable selection algorithms to identify, among
#' the variables in \code{xdata}, relevant predictors of the outcome variables
#' in \code{ydata}.
#'
#' @inheritParams SimulateGraphical
#' @param pk vector with the number of predictors in each independent block of
#'   variables in \code{xdata}. The number of independent blocks, which
#'   determines the maximum number of orthogonal latent variables that can be
#'   simulated, is given by \code{length(pk)}.
#' @param family type of outcome. If \code{family="gaussian"}, normally
#'   distributed outcomes are simulated. If \code{family="binomial"} or
#'   \code{family="multinomial"}, binary outcome(s) are simulated from a
#'   multinomial distribution where the probability is defined from a linear
#'   combination of normally distributed outcomes.
#' @param N number of classes of the categorical outcome. Only used if
#'   \code{family="multinomial"}.
#' @param ev_xz vector of the expected proportions of explained variances for each
#'   of the orthogonal latent variables. It must contain values in ]0,1[, and
#'   must be a vector of length \code{length(pk)} or a single value to generate
#'   latent variables with the same expected proportion of explained variance.
#' @param adjacency_x optional matrix encoding the conditional independence
#'   structure between predictor variables in \code{xdata}. This argument must
#'   be a binary symmetric matrix of size \code{sum(pk)} with zeros on the
#'   diagonal.
#' @param nu_within expected density (number of edges over the number of node
#'   pairs) of the conditional independence graph in the within-group blocks for
#'   predictors. For independent predictors, use \code{nu_within=0}. This
#'   argument is only used if \code{adjancency_x} is not provided.
#' @param theta_xz optional binary matrix encoding the predictor variables from
#'   \code{xdata} (columns) contributing to the definition of the orthogonal
#'   latent outcomes from \code{zdata} (rows).
#' @param nu_xz expected proportion of relevant predictors over the total number
#'   of predictors to be used for the simulation of the orthogonal latent
#'   outcomes. This argument is only used if \code{theta_xz} is not provided.
#' @param theta_zy optional binary matrix encoding the latent variables from
#'   \code{zdata} (columns) contributing to the definition of the observed
#'   outcomes from \code{ydata} (rows). This argument must be a square matrix of
#'   size \code{length(pk)}. If \code{theta_zy} is a diagonal matrix, each
#'   latent variable contributes to the definition of one observed outcome so
#'   that there is a one-to-one relationship between latent and observed
#'   outcomes (i.e. they are colinear). Nonzero off-diagonal elements in
#'   \code{theta_zy} introduce some correlation between the observed outcomes by
#'   construction from linear combinations implicating common latent outcomes.
#'   This argument is only used if \code{eta} is not provided.
#' @param nu_zy probability for each of the off-diagonal elements in
#'   \code{theta_zy} to be a 1. If \code{nu_zy=0}, \code{theta_zy} is a diagonal
#'   matrix. This argument is only used if \code{theta_zy} is not provided.
#' @param eta optional matrix of coefficients used in the linear combination of
#'   latent outcomes to generate observed outcomes.
#' @param eta_set vector defining the range of values from which \code{eta} is
#'   sampled. This argument is only used if \code{eta} is not provided.
#' @param ev_xx expected proportion of explained variance by the first Principal
#'   Component (PC1) of a Principal Component Analysis. This is the largest
#'   eigenvalue of the correlation (if \code{scale=TRUE}) or covariance (if
#'   \code{scale=FALSE}) matrix divided by the sum of eigenvalues. If
#'   \code{ev=NULL} (the default), the constant u is chosen by maximising the
#'   contrast of the correlation matrix.
#'
#' @return A list with: \item{xdata}{simulated predictor data.}
#'   \item{ydata}{simulated outcome data.} \item{proba}{simulated probability of
#'   belonging to each outcome class. Only used for \code{family="binomial"} or
#'   \code{family="multinomial"}.} \item{logit_proba}{logit of the simulated
#'   probability of belonging to each outcome class. Only used for
#'   \code{family="binomial"} or \code{family="multinomial"}.}
#'   \item{zdata}{simulated data for orthogonal latent outcomes.}
#'   \item{beta}{matrix of true beta coefficients used to generate outcomes in
#'   \code{ydata} from predictors in \code{xdata}.} \item{theta}{binary matrix
#'   indicating the predictors from \code{xdata} contributing to the definition
#'   of each of the outcome variables in \code{ydata}.} \item{eta}{matrix of
#'   coefficients used in the linear combination of latent variables from
#'   \code{zdata} to define observed outcomes in \code{ydata}.}
#'   \item{theta_zy}{binary matrix indicating the latent variables from
#'   \code{zdata} used in the definition of observed outcomes in \code{ydata}.}
#'   \item{xi}{matrix of true beta coefficients used to generate orthogonal
#'   latent outcomes in \code{zdata} from predictors in \code{xdata}.}
#'   \item{theta_xz}{binary matrix indicating the predictors from \code{xdata}
#'   contributing to the definition of each of the latent outcome variables in
#'   \code{zdata}.} \item{omega_xz}{precision matrix for variables in
#'   \code{xdata} and \code{zdata}.} \item{adjacency}{binary matrix encoding the
#'   conditional independence structure between variables from \code{xdata}
#'   (var), \code{zdata} (latent) and \code{ydata} (outcome).}
#'
#' @family simulation functions
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' ## Continuous outcomes
#'
#' # Univariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(5, 7, 3))
#' plot(Graph(simul$adjacency,
#'   satellites = TRUE,
#'   node_colour = c(rep("red", 3), rep("orange", 3), rep("skyblue", 15))
#' ))
#' print(simul$theta)
#'
#' # Multivariate outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(5, 7, 3))
#' plot(Graph(simul$adjacency,
#'   satellites = TRUE,
#'   node_colour = c(rep("red", 3), rep("orange", 3), rep("skyblue", 15))
#' ))
#' print(simul$theta)
#'
#' # Independent predictors
#' set.seed(1)
#' simul <- SimulateRegression(pk = c(5, 3), nu_within = 0)
#' plot(Graph(simul$adjacency,
#'   satellites = TRUE,
#'   node_colour = c(rep("red", 2), rep("orange", 2), rep("skyblue", 8))
#' ))
#'
#' # Blocks of strongly inter-connected predictors
#' set.seed(1)
#' simul <- SimulateRegression(
#'   pk = c(5, 5), nu_within = 0.5,
#'   v_within = c(0.5, 1), v_sign = -1, continuous = TRUE, pd_strategy = "min_eigenvalue"
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(simul$xdata),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' par(mar = rep(0, 4))
#' plot(Graph(simul$adjacency,
#'   satellites = TRUE,
#'   node_colour = c(rep("red", 2), rep("orange", 2), rep("skyblue", 10))
#' ))
#'
#'
#' ## Categorical outcomes
#'
#' # Binary outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 20, family = "binomial")
#' table(simul$ydata[, 1])
#'
#' # Categorical outcome
#' set.seed(1)
#' simul <- SimulateRegression(pk = 20, family = "multinomial")
#' apply(simul$ydata, 2, sum)
#' }
#' @export
SimulateRegression <- function(n = 100, pk = 10, N = 3,
                               family = "gaussian", ev_xz = 0.8,
                               adjacency_x = NULL, nu_within = 0.1,
                               theta_xz = NULL, nu_xz = 0.2,
                               theta_zy = NULL, nu_zy = 0.5,
                               eta = NULL, eta_set = c(-1, 1),
                               v_within = c(0.5, 1), v_sign = c(-1, 1), continuous = TRUE,
                               pd_strategy = "diagonally_dominant", ev_xx = NULL, scale = TRUE,
                               u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25) {
  # Checking the inputs
  if ((length(pk) > 1) & (family == "multinomial")) {
    stop("The simulation of multiple categorical outcomes is not possible with the current implementation.")
  }

  # Definition of the number of latent outcome variables
  q <- length(pk)
  p <- sum(pk)
  if (length(nu_xz) != q) {
    nu_xz <- rep(nu_xz[1], q)
  }
  if (length(ev_xz) != q) {
    ev_xz <- rep(ev_xz[1], q)
  }

  # Checking the values of ev_xz
  if (any(ev_xz <= 0) | any(ev_xz >= 1)) {
    stop("Invalid input for argument 'ev_xz'. Please provide values strictly between 0 and 1.")
  }

  # Simulation of the conditional independence structure with independent blocks
  if (is.null(adjacency_x)) {
    adjacency_x <- SimulateAdjacency(
      pk = pk, nu_between = 0, nu_within = nu_within,
      implementation = HugeAdjacency,
      topology = "random"
    )
  }

  # Simulation of the binary contribution status of predictors for latent outcome variables
  if (is.null(theta_xz)) {
    theta_xz <- SamplePredictors(pk = pk, q = q, nu = nu_xz, orthogonal = TRUE)
  }
  colnames(theta_xz) <- paste0("latent", 1:ncol(theta_xz))
  rownames(theta_xz) <- paste0("var", 1:nrow(theta_xz))

  # Simulation of precision matrix for both predictors and latent outcomes
  big_theta <- cbind(
    rbind(matrix(0, nrow = q, ncol = q), theta_xz),
    rbind(t(theta_xz), adjacency_x)
  )
  rownames(big_theta) <- colnames(big_theta)
  out <- SimulatePrecision(
    theta = big_theta, v_within = v_within,
    v_sign = v_sign, continuous = continuous,
    pd_strategy = pd_strategy, ev = ev_xx, scale = scale, u_list = u_list, tol = tol
  )
  omega <- out$omega

  # Setting diagonal precision for latent outcomes to reach expected proportion of explained variance
  xi <- NULL
  for (j in 1:q) {
    pred_ids <- seq(q + 1, q + p)
    omega[j, j] <- omega[j, pred_ids, drop = FALSE] %*% solve(omega[pred_ids, pred_ids]) %*% t(omega[j, pred_ids, drop = FALSE]) * 1 / ev_xz[j]
    xi <- cbind(xi, 1 / omega[j, j] * omega[j, pred_ids])
  }
  colnames(xi) <- colnames(theta_xz)

  # Computing the covariance matrix
  sigma <- solve(omega)

  # Simulation of data from multivariate normal distribution
  x <- MASS::mvrnorm(n, rep(0, p + q), sigma)
  rownames(x) <- paste0("obs", 1:nrow(x))

  # Separating predictors from latent outcome variables
  xdata <- x[, grep("var", colnames(x)), drop = FALSE]
  zdata <- x[, grep("latent", colnames(x)), drop = FALSE]

  if (family == "multinomial") {
    if (is.null(theta_zy)) {
      theta_zy <- SamplePredictors(pk = q, q = q, nu = nu_zy, orthogonal = FALSE)
    }

    proba <- matrix(NA, nrow = n, ncol = N)
    for (j in 1:N) {
      eta <- matrix(stats::runif(q * q, min = min(eta_set), max = max(eta_set)),
        ncol = q, nrow = q
      )
      eta <- eta * theta_zy

      rownames(eta) <- rownames(theta_zy) <- paste0("latent", 1:q)
      colnames(eta) <- colnames(theta_zy) <- paste0("outcome", 1:q)
      ydata <- zdata %*% eta

      if (j == 1) {
        beta <- xi %*% eta
        theta <- ifelse(beta != 0, yes = 1, no = 0)
      }

      proba[, j] <- 1 / (1 + exp(-ydata[, 1])) # inverse logit
    }
    proba <- t(apply(proba, 1, cumsum)) / N

    ydata_cat <- matrix(0, nrow = n, ncol = N)
    for (i in 1:n) {
      r <- stats::runif(n = 1)
      if (any(r < proba[i, ])) {
        ydata_cat[i, min(which(r < proba[i, ]))] <- 1
      }
    }

    # Setting row and column names
    rownames(ydata_cat) <- rownames(proba) <- rownames(ydata)
    colnames(ydata_cat) <- colnames(proba) <- paste0("cat", 1:N)
  } else {
    # Simulation of eta coefficients to get observed outcomes from latent outcomes
    if (is.null(eta)) {
      if (is.null(theta_zy)) {
        theta_zy <- SamplePredictors(pk = q, q = q, nu = nu_zy, orthogonal = FALSE)
      }
      eta <- matrix(stats::runif(q * q, min = min(eta_set), max = max(eta_set)),
        ncol = q, nrow = q
      )
      eta <- eta * theta_zy
    } else {
      theta_zy <- ifelse(eta != 0, yes = 1, no = 0)
    }
    rownames(eta) <- rownames(theta_zy) <- paste0("latent", 1:q)
    colnames(eta) <- colnames(theta_zy) <- paste0("outcome", 1:q)
    ydata <- zdata %*% eta

    # Computing the xy coefficients and binary contribution status
    beta <- xi %*% eta
    theta <- ifelse(beta != 0, yes = 1, no = 0)

    # Compute binary outcome for logistic regression
    if (family == "binomial") {
      proba <- matrix(NA, nrow = n, ncol = q)
      for (j in 1:q) {
        proba[, j] <- 1 / (1 + exp(-ydata[, j])) # inverse logit
      }

      ydata_cat <- matrix(0, nrow = n, ncol = q)
      for (j in 1:q) {
        for (i in 1:n) {
          ydata_cat[i, j] <- stats::rbinom(n = 1, size = 1, prob = proba[i, j])
        }
      }

      # Setting row and column names
      rownames(ydata_cat) <- rownames(proba) <- rownames(ydata)
      colnames(ydata_cat) <- colnames(proba) <- colnames(ydata)
    }
  }

  # Extracting the conditional independence structure between x, z and y
  adjacency <- rbind(
    cbind(matrix(0, nrow = q, ncol = q), t(theta_zy), matrix(0, nrow = q, ncol = p)),
    cbind(rbind(theta_zy, matrix(0, nrow = p, ncol = q)), big_theta)
  )
  rownames(adjacency) <- colnames(adjacency) <- c(colnames(theta_zy), rownames(big_theta))

  # Return the simulated X and Y
  if (family %in% c("binomial", "multinomial")) {
    if (family == "binomial") {
      out <- list(
        xdata = xdata, ydata = ydata_cat,
        proba = proba, logit_proba = ydata,
        zdata = zdata,
        beta = beta, theta = theta,
        eta = eta, theta_zy = theta_zy,
        xi = xi, theta_xz = theta_xz,
        omega_xz = omega,
        adjacency = adjacency
      )
    }
    if (family == "multinomial") {
      out <- list(
        xdata = xdata, ydata = ydata_cat,
        proba = proba, logit_proba = ydata,
        zdata = zdata,
        theta = theta,
        theta_zy = theta_zy,
        xi = xi, theta_xz = theta_xz,
        omega_xz = omega,
        adjacency = adjacency
      )
    }
  } else {
    out <- list(
      xdata = xdata, ydata = ydata, zdata = zdata,
      beta = beta, theta = theta,
      eta = eta, theta_zy = theta_zy,
      xi = xi, theta_xz = theta_xz,
      omega_xz = omega,
      adjacency = adjacency
    )
  }

  return(out)
}


#' Simulation of an undirected graph
#'
#' Simulates the adjacency matrix encoding an unweighted, undirected graph with
#' no self-loops.
#'
#' @inheritParams SimulateGraphical
#' @param pk number of nodes.
#' @param nu expected density (number of edges over the number of node pairs) of
#'   the graph. This argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}).
#' @param ... additional arguments to be passed to
#'   \code{\link[huge]{huge.generator}}.
#'
#' @return A symmetric adjacency matrix encoding an unweighted, undirected graph
#'   with no self-loops.
#'
#' @family simulation functions
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of a scale-free graph with 20 nodes
#' adjacency <- HugeAdjacency(pk = 20, topology = "scale-free")
#' plot(Graph(adjacency))
#' }
#' @export
HugeAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
  # Storing extra arguments
  extra_args <- list(...)

  # Extracting relevant extra arguments
  ids <- which(names(extra_args) %in% names(formals(huge::huge.generator)))
  ids <- ids[!ids %in% c("n", "d", "prob", "graph", "verbose")]

  # Running simulation model
  mymodel <- do.call(huge::huge.generator, args = c(
    list(
      n = 2, d = sum(pk), prob = nu,
      graph = topology, verbose = FALSE
    ),
    extra_args[ids]
  ))
  theta <- as.matrix(mymodel$theta)

  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids <- sample(ncol(theta))
  theta <- theta[ids, ids]

  return(theta)
}


#' Simulation of an undirected graph with block structure
#'
#' Simulates the adjacency matrix of an unweighted, undirected graph with no
#' self-loops, and with different densities in diagonal compared to off-diagonal
#' blocks.
#'
#' @inheritParams SimulateGraphical
#'
#' @return A symmetric adjacency matrix encoding an unweighted, undirected graph
#'   with no self-loops, and with different densities in diagonal compared to off-diagonal
#' blocks.
#'
#' @family simulation functions
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of a scale-free graph with 20 nodes
#' adjacency <- SimulateAdjacency(pk = 20, topology = "scale-free")
#' plot(Graph(adjacency))
#'
#' # Simulation of a random graph with block structure
#' adjacency <- SimulateAdjacency(
#'   pk = rep(10, 3),
#'   nu_within = 0.7, nu_between = 0.03
#' )
#' plot(Graph(adjacency))
#'
#' # User-defined function for graph simulation
#' CentralNode <- function(pk, hub = 1) {
#'   theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
#'   theta[hub, ] <- 1
#'   theta[, hub] <- 1
#'   diag(theta) <- 0
#'   return(theta)
#' }
#' simul <- SimulateAdjacency(pk = 10, implementation = CentralNode)
#' plot(Graph(simul)) # star
#' simul <- SimulateAdjacency(pk = 10, implementation = CentralNode, hub = 2)
#' plot(Graph(simul)) # variable 2 is the central node
#' }
#'
#' @export
SimulateAdjacency <- function(pk = 10,
                              implementation = HugeAdjacency, topology = "random",
                              nu_within = 0.1, nu_between = 0, ...) {
  # Storing all arguments
  args <- c(mget(ls()), list(...))

  # Checking the inputs
  if (topology != "random") {
    if (length(pk) > 1) {
      pk <- sum(pk)
      warning(paste0("Multi-block simulations are only allowed with topology='random'. Argument 'pk' has been set to ", pk, "."))
    }
  }

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Identifying relevant arguments
  if (!"..." %in% names(formals(implementation))) {
    ids <- which(names(args) %in% names(formals(implementation)))
    args <- args[ids]
  }

  # Simulation of the adjacency matrix
  if ("nu" %in% names(formals(implementation))) {
    if (length(pk) > 1) {
      # Initialising theta
      theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
      theta_vect <- theta[upper.tri(theta)]

      # Allowing for different densities in within and between blocks
      theta_w <- do.call(implementation, args = c(args, list(nu = nu_within)))
      theta_w_vect <- theta_w[upper.tri(theta_w)]
      theta_b <- do.call(implementation, args = c(args, list(nu = nu_between)))
      theta_b_vect <- theta_b[upper.tri(theta_b)]

      # Filling within and between blocks
      for (k in block_ids) {
        if (k %in% unique(diag(bigblocks))) {
          theta_vect[bigblocks_vect == k] <- theta_w_vect[bigblocks_vect == k]
        } else {
          theta_vect[bigblocks_vect == k] <- theta_b_vect[bigblocks_vect == k]
        }
      }
      theta[upper.tri(theta)] <- theta_vect
      theta <- theta + t(theta)
    } else {
      theta <- do.call(implementation, args = c(args, list(nu = nu_within)))
    }
  } else {
    theta <- do.call(implementation, args = c(args))
  }

  # Ensuring the adjacency matrix is symmetric (undirected graph) with no self-loops
  theta <- ifelse(theta + t(theta) != 0, yes = 1, no = 0)
  diag(theta) <- 0

  # Setting variable names
  colnames(theta) <- rownames(theta) <- paste0("var", 1:ncol(theta))

  return(theta)
}


#' Simulation of symmetric matrix with block structure
#'
#' Simulates a symmetric matrix with block structure. Matrix entries are
#' sampled from (i) a discrete uniform distribution taking values in
#' \code{v_within} (for entries in the diagonal block) or \code{v_between} (for
#' entries in off-diagonal blocks) if \code{continuous=FALSE}, or (ii) a
#' continuous uniform distribution taking values in the range given by
#' \code{v_within} or \code{v_between} if \code{continuous=TRUE}.
#'
#' @param pk vector of the number of variables per group, defining the block
#'   structure.
#' @param v_within vector defining the (range of) nonzero entries in the
#'   diagonal blocks. If \code{continuous=FALSE}, \code{v_within} is the set of
#'   possible values. If \code{continuous=FALSE}, \code{v_within} is the range
#'   of possible values.
#' @param v_between vector defining the (range of) nonzero entries in the
#'   off-diagonal blocks. If \code{continuous=FALSE}, \code{v_between} is the
#'   set of possible precision values. If \code{continuous=FALSE},
#'   \code{v_between} is the range of possible precision values. This argument
#'   is only used if \code{length(pk)>1}.
#' @param v_sign vector of possible signs for matrix entries. Possible
#'   inputs are: \code{-1} for negative entries only, \code{1} for
#'   positive entries only, or \code{c(-1, 1)} for both positive and
#'   negative entries.
#' @param continuous logical indicating whether to sample precision values from
#'   a uniform distribution between the minimum and maximum values in
#'   \code{v_within} (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=TRUE}) or from proposed values in \code{v_within}
#'   (diagonal blocks) or \code{v_between} (off-diagonal blocks)
#'   (\code{continuous=FALSE}).
#'
#' @return A symmetric matrix with uniformly distributed entries sampled from
#'   different distributions for diagonal and off-diagonal blocks.
#'
#' @export
SimulateSymmetricMatrix <- function(pk = 10,
                                    v_within = c(0.5, 1), v_between = c(0, 0.1),
                                    v_sign = c(-1, 1), continuous = FALSE) {
  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Building absolute v matrix
  v <- bigblocks
  v_vect <- v[upper.tri(v)]
  for (k in block_ids) {
    if (k %in% v_vect) {
      if (k %in% unique(diag(bigblocks))) {
        if (continuous) {
          v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_within), max = max(v_within))
        } else {
          v_vect[bigblocks_vect == k] <- base::sample(v_within, size = sum(bigblocks_vect == k), replace = TRUE)
        }
      } else {
        if (continuous) {
          v_vect[bigblocks_vect == k] <- stats::runif(sum(bigblocks_vect == k), min = min(v_between), max = max(v_between))
        } else {
          v_vect[bigblocks_vect == k] <- base::sample(v_between, size = sum(bigblocks_vect == k), replace = TRUE)
        }
      }
    }
  }

  # Sampling the sign of precision entries
  v_vect <- v_vect * base::sample(sort(unique(v_sign)), size = length(v_vect), replace = TRUE)

  # Building v matrix
  diag(v) <- 0
  v[upper.tri(v)] <- v_vect
  v[lower.tri(v)] <- 0
  v <- v + t(v)

  return(v)
}


#' Simulation of a precision matrix
#'
#' Simulates a sparse precision matrix from an adjacency matrix \code{theta}
#' encoding a conditional independence graph. Zero entries in the precision
#' matrix indicate pairwise conditional independence. Diagonal entries can be
#' tuned to (i) maximise the contrast of the correlation matrix, or (ii) reach a
#' user-defined proportion of explained variance by the first Principal
#' Component (see \code{\link{MakePositiveDefinite}}).
#'
#' @inheritParams SimulateGraphical
#' @param theta binary and symmetric adjacency matrix encoding the conditional
#'   independence structure.
#' @param scale logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix.
#'
#' @return A list with: \item{omega}{true simulated precision matrix.}
#'   \item{u}{value of the constant u used to ensure that \code{omega} is
#'   positive definite.}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of an adjacency matrix
#' theta <- SimulateAdjacency(pk = c(5, 5), nu_within = 0.7)
#' print(theta)
#'
#' # Simulation of a precision matrix maximising the contrast
#' simul <- SimulatePrecision(theta = theta)
#' print(simul$omega)
#'
#' # Simulation of a precision matrix with specific ev by PC1
#' simul <- SimulatePrecision(
#'   theta = theta,
#'   pd_strategy = "min_eigenvalue",
#'   ev = 0.3, scale = TRUE
#' )
#' print(simul$omega)
#' }
#' @export
SimulatePrecision <- function(pk = NULL, theta,
                              v_within = c(0.5, 1), v_between = c(0, 0.1),
                              v_sign = c(-1, 1), continuous = TRUE,
                              pd_strategy = "diagonally_dominant", ev = NULL, scale = TRUE,
                              u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25) {
  # Checking inputs and defining pk
  if (is.null(pk)) {
    pk <- ncol(theta)
  } else {
    if (sum(pk) != ncol(theta)) {
      stop("Arguments 'pk' and 'theta' are not consistent. The sum of 'pk' entries must be equal to the number of rows and columns in 'theta'.")
    }
  }

  # Checking the choice of pd_strategy
  if (!pd_strategy %in% c("diagonally_dominant", "min_eigenvalue")) {
    stop("Invalid input for argument 'pd_strategy'. Possible values are: 'diagonally_dominant' or 'min_eigenvalue'.")
  }

  # Checking other input values
  if (any((v_within < 0) | (v_within > 1))) {
    stop("Invalid input for argument 'v_within'. Values must be between 0 and 1.")
  }
  if (any((v_between < 0) | (v_between > 1))) {
    stop("Invalid input for argument 'v_between'. Values must be between 0 and 1.")
  }
  if (any(!v_sign %in% c(-1, 1))) {
    stop("Invalid input for argument 'v_sign'. Possible values are -1 and 1.")
  }

  # Ensuring that v values are lower than or equal to 1
  if (any(abs(v_within) > 1)) {
    v_within <- v_within / max(abs(v_within))
    message("The values provided in 'v_within' have been re-scaled to be lower than or equal to 1 in absolute value.")
  }

  # Ensuring that diagonal entries of theta are zero
  diag(theta) <- 0

  # Building v matrix
  v <- SimulateSymmetricMatrix(
    pk = pk, v_within = v_within, v_between = v_between,
    v_sign = v_sign, continuous = continuous
  )

  # Filling off-diagonal entries of the precision matrix
  omega_tilde <- theta * v

  # Ensuring positive definiteness
  omega_pd <- MakePositiveDefinite(
    omega = omega_tilde, pd_strategy = pd_strategy,
    ev = ev, scale = scale, u_list = u_list, tol = tol
  )

  # # Preparing realistic diagonally dominant precision matrix
  # if (pd_strategy == "diagonally_dominant") {
  #   # Defining grid of u values if not provided
  #   if (is.null(u)) {
  #     u <- 10^-(seq(0, 5, by = 0.1))
  #     refining_u_grid <- TRUE
  #     niter_max <- 5
  #     tolerance <- 10
  #   } else {
  #     refining_u_grid <- FALSE
  #   }
  #
  #   # Filling off-diagonal entries of the precision matrix
  #   omega <- theta * v
  #
  #   # Calibrate u based on contrasts of the correlation matrix
  #   contrast <- NULL
  #   for (u_value in u) {
  #     omega_tmp <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
  #     C <- stats::cov2cor(solve(omega_tmp))
  #     contrast <- c(contrast, Contrast(C))
  #   }
  #
  #   # Avoiding extreme values in u grid if not provided by the user
  #   if (refining_u_grid) {
  #     stop <- 0
  #     niter <- 1
  #     while (stop == 0) {
  #       niter <- niter + 1
  #       if (niter == niter_max_u_grid) {
  #         stop <- 1
  #       }
  #       # Satisfied with calibrated u if the argmax is not too close to the boundaries (as defined from tolerance_u_grid)
  #       if (any(which(contrast == max(contrast)) %in% seq(tolerance_u_grid, length(u) - tolerance_u_grid) == TRUE)) {
  #         stop <- 1
  #       } else {
  #         # Adding smaller values of u
  #         if (any(which(contrast == max(contrast)) %in% seq(1, tolerance_u_grid) == TRUE)) {
  #           u <- c(u, 10^-seq(min(-log10(u)) - u_delta, min(-log10(u)), by = 0.1))
  #         }
  #
  #         # Adding larger values of u
  #         if (any(which(contrast == max(contrast)) %in% seq(length(u) - tolerance_u_grid, length(u)) == TRUE)) {
  #           u <- c(u, 10^-seq(max(-log10(u)), max(-log10(u) + u_delta), by = 0.1))
  #         }
  #
  #         # Sorting values in u
  #         u <- sort(u, decreasing = TRUE)
  #
  #         # Computing the contrast for all visited values of u
  #         contrast <- NULL
  #         for (u_value in u) {
  #           omega_tmp <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
  #           C <- stats::cov2cor(solve(omega_tmp))
  #           contrast <- c(contrast, Contrast(C))
  #         }
  #       }
  #     }
  #   }
  #
  #   # Computing calibrated precision matrix
  #   if (length(u) > 1) {
  #     u_value <- u[length(contrast) - which.max(rev(contrast)) + 1] # adding smallest possible u value to the diagonal
  #     omega <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
  #   } else {
  #     omega <- omega_tmp
  #   }
  # }
  #
  # # Allowing for higher correlations using smallest eigenvalue
  # if (pd_strategy == "min_eigenvalue") {
  #   # Defining a small constant
  #   if (is.null(u)) {
  #     u <- 1e-5
  #   }
  #
  #   # Initialisation of full precision matrix
  #   omega <- matrix(0, ncol = sum(pk), nrow = sum(pk))
  #
  #   # Creating matrix with block indices
  #   bigblocks <- BlockMatrix(pk)
  #
  #   # Making positive definite diagonal blocks (Schur complement shows it is p.d.)
  #   for (i in unique(diag(bigblocks))) {
  #     # Filling off-diagonal entries of the precision matrix (for corresponding diagonal block)
  #     omega_block <- matrix((theta * v)[which(bigblocks == i)], ncol = sum(diag(bigblocks) == i))
  #
  #     # # Computing the signed adjacency
  #     # signed_theta <- sign(omega_block)
  #
  #     # Defining a positive definite precision matrix
  #     omega_tmp <- MakePositiveDefinite(omega = omega_block, u_value = u, pd_strategy = pd_strategy)
  #
  #     # # Using sampled entries (need to be lower than or equal to 1 in absolute value)
  #     # diag(omega_block) <- diag(omega_tmp)
  #
  #     # Filling full precision matrix
  #     omega[which(bigblocks == i)] <- omega_tmp
  #   }
  #
  #   # Accounting for off-diagonal blocks (sum of p.d. is p.d.)
  #   if (length(pk) > 1) {
  #     # Filling off-diagonal entries of the precision matrix (for off-diagonal blocks)
  #     omega_block <- theta * v
  #     omega_block[which(bigblocks %in% unique(diag(bigblocks)))] <- 0
  #
  #     # # Computing the signed adjacency
  #     # signed_theta <- sign(omega_block)
  #
  #     # Defining a positive definite precision matrix
  #     omega_tmp <- MakePositiveDefinite(omega = omega_block, u_value = u, pd_strategy = pd_strategy)
  #
  #     # # Using sampled entries (need to be lower than or equal to 1 in absolute value)
  #     # diag(omega_block) <- diag(omega_tmp)
  #
  #     # Computing sum of the p.d. matrices
  #     omega <- omega + omega_tmp
  #   }
  #
  #   # Setting row and column names
  #   rownames(omega) <- colnames(omega) <- colnames(theta)
  # }

  # Returning the output
  return(omega_pd)
}


#' Simulation of binary contribution status
#'
#' Simulates the binary contribution status of potential predictor variables
#' from different blocks to outcome variables. For each outcome, the set of true
#' predictors is sampled from one block of potential predictors. If the blocks
#' of variables are independent, the outcomes will be independent too.
#'
#' @inheritParams SimulateSymmetricMatrix
#' @param q number of outcome variables. By default, one block of predictor is
#'   linked to one outcome, i.e. \code{q=sum(pk)}.
#' @param nu vector of probabilities. Each entry corresponds to one block of
#'   predictors and defines the probability for each predictor within the block
#'   to be chosen as true predictor of the corresponding outcome variable.
#' @param orthogonal logical indicating if the outcomes have to be defined from
#'   independent blocks of predictors as encoded in \code{pk}.
#'
#' @return A binary matrix encoding the contribution status of each predictor
#'   variable (columns) to each outcome variable (rows).
#'
#' @export
SamplePredictors <- function(pk, q = NULL, nu = 0.1, orthogonal = TRUE) {
  # Definition of the number of outcome variables
  if (is.null(q)) {
    q <- length(pk)
  }
  if (length(nu) != q) {
    nu <- rep(nu[1], q)
  }

  # Simulation of the binary status for true predictors
  theta <- matrix(0, nrow = q, ncol = sum(pk))
  for (k in 1:q) {
    if (orthogonal) {
      if (k > 1) {
        ids <- seq(cumsum(pk)[k - 1] + 1, cumsum(pk)[k])
      } else {
        ids <- seq(1, cumsum(pk)[k])
      }
      theta[k, ids] <- stats::rbinom(pk[k], size = 1, prob = nu[k])

      # Introducing at least one true predictor
      if (sum(theta[k, ids]) == 0) {
        theta[k, sample(ids, size = 1)] <- 1
      }
    } else {
      theta[k, ] <- stats::rbinom(sum(pk), size = 1, prob = nu[k])
      theta[k, k] <- 1
    }
  }

  return(t(theta))
}


#' Making positive definite
#'
#' Determines the diagonal entries of a symmetric matrix to ensure it is
#' positive definite. For this, diagonal entries of the matrix are defined to be
#' higher than (i) the sum of entries on the corresponding rows, which ensure it
#' is diagonally dominant, or (ii) the absolute value of the smallest eigenvalue
#' of the same matrix with a diagonal of zeros. The magnitude of (standardised)
#' values in the inverse matrix is tuned by adding a constant u to the diagonal
#' entries. Considering the matrix to make positive definite is a precision
#' matrix, the constant u is chosen to (i) maximise the contrast of the
#' corresponding correlation matrix, or (ii) tune the proportion of explained
#' variance by the first Principal Component (i.e. largest eigenvalue of the
#' covariance/correlation matrix divided by the sum of eigenvalues).
#'
#' @inheritParams SimulateGraphical
#' @param omega input matrix.
#' @param scale logical indicating if the proportion of explained variance by
#'   PC1 should be computed from the correlation (\code{scale=TRUE}) or
#'   covariance (\code{scale=FALSE}) matrix.
#'
#' @return A list with: \item{omega}{positive definite matrix.} \item{u}{value
#'   of the constant u.}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of a symmetric matrix
#' p <- 5
#' set.seed(1)
#' omega <- matrix(rnorm(p * p), ncol = p)
#' omega <- omega + t(omega)
#' diag(omega) <- 0
#'
#' # Diagonal dominance maximising contrast
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant"
#' )
#' eigen(omega_pd$omega)$values # positive eigenvalues
#'
#' # Diagonal dominance with specific proportion of explained variance by PC1
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant",
#'   ev = 0.55
#' )
#' lambda_inv <- 1 / eigen(omega_pd$omega)$values
#' max(lambda_inv) / sum(lambda_inv) # expected ev
#'
#' # Diagonal dominance with specific proportion of explained variance by PC1 (scaled)
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "diagonally_dominant",
#'   ev = 0.55, scale = TRUE
#' )
#' lambda_inv <- eigen(cov2cor(solve(omega_pd$omega)))$values
#' max(lambda_inv) / sum(lambda_inv) # expected ev
#'
#' # Non-negative eigenvalues maximising contrast
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue"
#' )
#' eigen(omega_pd$omega)$values # positive eigenvalues
#'
#' # Non-negative eigenvalues with specific proportion of explained variance by PC1
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue",
#'   ev = 0.7
#' )
#' lambda_inv <- 1 / eigen(omega_pd$omega)$values
#' max(lambda_inv) / sum(lambda_inv)
#'
#' # Non-negative eigenvalues with specific proportion of explained variance by PC1 (scaled)
#' omega_pd <- MakePositiveDefinite(omega,
#'   pd_strategy = "min_eigenvalue",
#'   ev = 0.7, scale = TRUE
#' )
#' lambda_inv <- eigen(cov2cor(solve(omega_pd$omega)))$values
#' max(lambda_inv) / sum(lambda_inv)
#' }
#' @export
MakePositiveDefinite <- function(omega, pd_strategy = "diagonally_dominant",
                                 ev = NULL, scale = FALSE, u_list = c(1e-10, 1),
                                 tol = .Machine$double.eps^0.25) {

  # Making positive definite using diagonally dominance
  if (pd_strategy == "diagonally_dominant") {
    # Constructing the diagonal as the sum of entries
    diag(omega) <- apply(abs(omega), 1, sum)
    lambda <- eigen(omega)$values
  }
  # Making positive definite using eigendecomposition
  if (pd_strategy == "min_eigenvalue") {
    # Extracting smallest eigenvalue of omega_tilde
    lambda <- eigen(omega)$values
    lambda0 <- abs(min(lambda))

    # Making the precision matrix positive semidefinite
    lambda <- lambda + lambda0
    diag(omega) <- lambda0
  }

  if (is.null(ev)) {
    # Finding u that maximises the contrast
    argmax_u <- stats::optimise(MaxContrast,
      omega = omega, maximum = TRUE, tol = 1,
      lower = min(u_list), upper = max(u_list)
    )
    u <- argmax_u$maximum
  } else {
    # Finding extreme values
    if (scale) {
      max_ev <- TuneExplainedVarianceCor(u = min(u_list), omega = omega)
      min_ev <- TuneExplainedVarianceCor(u = max(u_list), omega = omega)
    } else {
      max_ev <- TuneExplainedVarianceCov(u = min(u_list), lambda = lambda)
      min_ev <- TuneExplainedVarianceCov(u = max(u_list), lambda = lambda)
    }

    # Finding u corresponding to the required proportion of explained variance
    if ((ev <= min_ev) | (ev >= max_ev)) {
      if (ev <= min_ev) {
        u <- max(u_list)
        if (ev < min_ev) {
          message(paste0("The smallest proportion of explained variance by PC1 that can be obtained is ", round(min_ev, digits = 2), "."))
        }
      } else {
        u <- min(u_list)
        if (ev > max_ev) {
          message(paste0("The largest proportion of explained variance by PC1 that can be obtained is ", round(max_ev, digits = 2), "."))
        }
      }
    } else {
      if (scale) {
        # Minimising the difference between requested and possible ev
        argmin_u <- stats::optimise(TuneExplainedVarianceCor,
          omega = omega, ev = ev,
          lower = min(u_list), upper = max(u_list), tol = tol
        )
        u <- argmin_u$minimum
      } else {
        # Minimising the difference between requested and possible ev
        argmin_u <- stats::optimise(TuneExplainedVarianceCov,
          lambda = lambda, ev = ev,
          lower = min(u_list), upper = max(u_list), tol = tol
        )
        u <- argmin_u$minimum
      }
    }
  }

  # Constructing the diagonal
  diag(omega) <- diag(omega) + u

  return(list(omega = omega, u = u))
}


#' Maximising matrix contrast
#'
#' Computes the contrast of the correlation matrix obtained by adding u to the
#' diagonal of the precision matrix. This function is used to find the value of
#' u that maximises the contrast when constructing a diagonally dominant
#' precision matrix.
#'
#' @param u constant u added to the diagonal of the precision matrix.
#' @param omega positive semi-definite precision matrix.
#' @param digits number of digits to use in the definition of the contrast.
#'
#' @return A single number, the contrast of the generated precision matrix.
MaxContrast <- function(u, omega, digits = 3) {
  diag(omega) <- diag(omega) + u
  return(Contrast(stats::cov2cor(omega), digits = digits))
}


#' Matrix contrast
#'
#' Computes the matrix contrast, defined as the number of unique truncated
#' entries with a specified number of digits.
#'
#' @param mat input matrix.
#' @param digits number of digits to use.
#'
#' @return A single number, the contrast of the input matrix.
#'
#' @export
Contrast <- function(mat, digits = 3) {
  return(length(unique(round(as.vector(abs(mat)), digits = digits))))
}


#' Tuning function (covariance)
#'
#' Computes the difference in absolute value between the desired and observed
#' proportion of explained variance from the first Principal Component of a
#' Principal Component Analysis applied on the covariance matrix. The precision
#' matrix is obtained by adding u to the diagonal of a positive semidefinite
#' matrix. This function is used to find the value of the constant u
#' that generates a covariance matrix with desired proportion of explained
#' variance.
#'
#' @inheritParams MaxContrast
#' @param ev desired proportion of explained variance. If \code{ev=NULL}, the
#'   obtained proportion of explained variance is returned.
#' @param lambda eigenvalues of the positive semidefinite precision matrix.
#'
#' @return The difference in proportion of explained variance in absolute values
#'   or observed proportion of explained variance (if \code{ev=NULL}).
TuneExplainedVarianceCov <- function(u, ev = NULL, lambda) {
  lambda <- lambda + u
  lambda_inv <- 1 / lambda
  tmp_ev <- max(lambda_inv) / sum(lambda_inv)
  if (is.null(ev)) {
    out <- tmp_ev
  } else {
    out <- abs(tmp_ev - ev)
  }
  return(out)
}


#' Tuning function (correlation)
#'
#' Computes the difference in absolute value between the desired and observed
#' proportion of explained variance from the first Principal Component of a
#' Principal Component Analysis applied on the correlation matrix. The precision
#' matrix is obtained by adding u to the diagonal of a positive semidefinite
#' matrix. This function is used to find the value of the constant u
#' that generates a correlation matrix with desired proportion of explained
#' variance.
#'
#' @inheritParams TuneExplainedVarianceCov
#' @param omega positive semidefinite precision matrix.
#'
#' @return The difference in proportion of explained variance in absolute values
#'   or observed proportion of explained variance (if \code{ev=NULL}).
TuneExplainedVarianceCor <- function(u, ev = NULL, omega) {
  diag(omega) <- diag(omega) + u
  mycor <- stats::cov2cor(solve(omega))
  tmp_ev <- norm(mycor, type = "2") / ncol(mycor)
  if (is.null(ev)) {
    out <- tmp_ev
  } else {
    out <- abs(tmp_ev - ev)
  }
  return(out)
}

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
#'   functionalities implemented in \code{\link[huge]{huge.generator}} are used.
#'   Alternatively, a user-defined function can be used. It must take \code{pk},
#'   \code{topology} and \code{nu} as arguments and return a
#'   \code{(sum(pk)*(sum(pk)))} binary and symmetric matrix for which diagonal
#'   entries are all equal to zero. This function is only applied if
#'   \code{theta} is not provided.
#' @param topology topology of the simulated graph. If using
#'   \code{implementation=SimulateAdjacency}, possible values are listed for the
#'   argument \code{graph} of \code{\link[huge]{huge.generator}}. These are:
#'   "random", "hub", "cluster", "band" and "scale-free".
#' @param nu_within expected density of within-group blocks the graph. If
#'   \code{length(pk)=1}, this is the expected density of the graph. If
#'   \code{implementation=SimulateAdjacency}, this argument is only used for
#'   \code{topology="random"} or \code{topology="cluster"} (see argument
#'   \code{prob} in \code{\link[huge]{huge.generator}}).
#' @param nu_between expected density of between-group blocks in the graph.
#'   Similar to \code{nu_within}. By default, the same density is used for
#'   within and between blocks (\code{nu_within}=\code{nu_between}). Only used
#'   if \code{length(pk)>1}.
#' @param output_matrices logical indicating if the true precision and (partial)
#'   correlation matrices should be included in the output.
#' @param v_within vector defining the (range of) nonzero entries in the
#'   diagonal blocks of the precision matrix. If \code{continuous=FALSE},
#'   \code{v_within} is the set of possible precision values. If
#'   \code{continuous=FALSE}, \code{v_within} is the range of possible precision
#'   values.
#' @param v_between vector defining the (range of) nonzero entries in the
#'   off-diagonal blocks of the precision matrix. If \code{continuous=FALSE},
#'   \code{v_between} is the set of possible precision values. If
#'   \code{continuous=FALSE}, \code{v_between} is the range of possible
#'   precision values. This argument is only used if \code{length(pk)>1}.
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
#'   \code{pd_strategy="nonnegative_eigenvalues"}, diagonal entries are set to
#'   the sum of the absolute value of the smallest eigenvalue and a constant u.
#' @param u optional vector of values for constant u used to ensure positive
#'   definiteness of the simulated precision matrix. The value that maximises
#'   the contrast of the simulated correlation matrix over the grid \code{u} is
#'   used. If \code{u=NULL}, a grid of values is automatically generated and
#'   iteratively shifted to ensure that the chosen value is not an extreme
#'   value.
#' @param niter_max_u_grid maximum number of iterations where the grid of u
#'   values is shifted. This parameter is only used with \code{u=NULL}.
#' @param tolerance_u_grid number of values between the chosen value for u and
#'   the end of the grid (first or last value). This parameter is only used with
#'   \code{u=NULL}.
#' @param u_delta difference in log10-scale between the smallest (or largest) u
#'   values in the grid used in the current iteration and the one built for next
#'   iteration. This parameter is only used with \code{u=NULL}.
#' @param ... additional arguments passed to the graph simulation function
#'   provided in \code{implementation}.
#'
#' @seealso \code{\link{SimulatePrecision}}, \code{\link{MakePositiveDefinite}},
#'   \code{\link{Contrast}}, \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables.} \item{theta}{adjacency matrix of the simulated
#'   graph} \item{omega}{true simulated precision matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{phi}{true simulated partial correlation
#'   matrix. Only returned if \code{output_matrices=TRUE}.} \item{C}{true
#'   simulated correlation matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u}{chosen value of u. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u_grid}{grid of u values. Only returned
#'   if \code{output_matrices=TRUE}.} \item{contrast_path}{contrast values
#'   obtained for the values of u listed in u_grid. Only returned if
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
#' simul <- SimulateGraphical(n = 100, pk = pk, nu_within = 0.05)
#' mycor <- cor(simul$data)
#' Heatmap(mycor,
#'   colours = c("darkblue", "white", "firebrick3"),
#'   legend_range = c(-1, 1), legend_length = 50, legend = FALSE
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
#' CentralNode <- function(pk, topology = NULL, nu = NULL, hub = 1) {
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
                              implementation = SimulateAdjacency, topology = "random",
                              nu_within = 0.1, nu_between = NULL,
                              output_matrices = FALSE,
                              v_within = c(-1, 1), v_between = c(-0.1, 0.1), continuous = FALSE,
                              pd_strategy = "diagonally_dominant",
                              u = NULL, niter_max_u_grid = 5, tolerance_u_grid = 10, u_delta = 5, ...) {
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
    theta <- SimulateBlockAdjacency(
      pk = pk,
      implementation = implementation, topology = topology,
      nu_within = nu_within, nu_between = nu_between, ...
    )
  }

  # Simulation of a precision matrix
  out <- SimulatePrecision(
    pk = pk, theta = theta,
    v_within = v_within, v_between = v_between, continuous = continuous,
    pd_strategy = pd_strategy, u = u,
    niter_max_u_grid = niter_max_u_grid, tolerance_u_grid = tolerance_u_grid, u_delta = u_delta
  )
  omega <- out$omega

  # Computing the correlation matrix
  C <- stats::cov2cor(solve(omega)) # true correlation matrix

  # Computing the partial correlation matrix
  if (output_matrices) {
    phi <- -stats::cov2cor(omega) + 2 * diag(ncol(omega))
  }

  # Simulating data from multivariate normal distribution
  x <- MASS::mvrnorm(n, rep(0, p), C)
  colnames(x) <- paste0("var", 1:ncol(x))
  rownames(x) <- paste0("obs", 1:nrow(x))

  if (output_matrices) {
    return(list(
      data = x, theta = theta,
      omega = omega, phi = phi, C = C,
      u = out$u, u_grid = out$u_grid, contrast_path = out$contrast_path
    ))
  } else {
    return(list(data = x, theta = theta))
  }
}


#' Simulation of data with underlying clusters
#'
#' Simulates multivariate Normal data with clusters of participants sharing
#' similar variable profiles. This simulator is based on a graph structure
#' encoding conditional correlations between the observations.
#'
#' @inheritParams SimulateGraphical
#' @param n vector of the number of observations per cluster in the simulated
#'   data. The number of observations in the simulated data is \code{sum(n)}.
#' @param pk vector of the number of variables in the simulated data.
#' @param adjacency optional binary and symmetric adjacency matrix encoding the
#'   conditional graph structure between observations. The clusters encoded in
#'   this argument must be in line with those indicated in \code{n}. To generate
#'   a block structure, the within-cluster observations must be more likely to
#'   be correlated than between clusters. For between-cluster relationships to
#'   be apparent, \code{v_between} must be nonzero. To be used with care.
#'
#' @seealso \code{\link{MakePositiveDefinite}}, \code{\link{Contrast}},
#'   \code{\link{GraphicalModel}}
#' @family simulation functions
#'
#' @return A list with: \item{data}{simulated data with \code{n} observation and
#'   \code{sum(pk)} variables} \item{theta}{adjacency matrix of the simulated
#'   graph} \item{omega}{true simulated precision matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{phi}{true simulated partial correlation
#'   matrix. Only returned if \code{output_matrices=TRUE}.} \item{C}{true
#'   simulated correlation matrix. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u}{chosen value of u. Only returned if
#'   \code{output_matrices=TRUE}.} \item{u_grid}{grid of u values. Only returned
#'   if \code{output_matrices=TRUE}.} \item{contrast_path}{contrast values
#'   obtained for the values of u listed in u_grid. Only returned if
#'   \code{output_matrices=TRUE}.}
#'
#' @details By default, with \code{nu=1}, a graph with \code{length(n)} fully
#'   connected components is used to infer the covariance matrix (between
#'   observations). The data is simulated from the multivariate Normal
#'   distribution with the generated covariance. This is based on similar
#'   computations to \code{\link{SimulateGraphical}} but applied on the
#'   transpose of the data matrix (i.e. encoding relationships between
#'   observations instead of variables). Using a smaller value for \code{nu} or
#'   a larger range for \code{v_within} will generate more heterogeneous blocks
#'   of correlated observations. Note that as \code{v_within} is controlling
#'   entries in the precision matrix, it must be negative to yield positive
#'   correlations.
#'
#' @family simulation functions
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of 15 observations belonging to 3 groups
#' set.seed(1)
#' simul <- SimulateClustering(n = c(5, 5, 5), pk = 100)
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # Simulation with weaker within-cluster correlations
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100,
#'   v_within = c(-1, -0.1), continuous = TRUE
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # Introducing between-cluster correlations
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100, nu_within = 1, nu_between = 0.3,
#'   v_within = c(-1, -0.5), v_between = c(-0.5, 0), continuous = TRUE
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # User-defined structure: between-cluster pair obs1-obs8
#' adjacency <- CoMembership(c(rep(1, 5), rep(2, 5), rep(3, 5)))
#' adjacency[1, 8] <- adjacency[8, 1] <- 1
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100, adjacency = adjacency, output_matrices = TRUE,
#'   v_within = c(-1, -0.5), v_between = -1, continuous = TRUE
#' )
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#' }
#' @export
SimulateClustering <- function(n = c(10, 10), pk = 20, adjacency = NULL,
                               nu_within = 1, nu_between = 1, output_matrices = FALSE,
                               v_within = c(-1, -0.5), v_between = 0, continuous = TRUE,
                               pd_strategy = "diagonally_dominant",
                               u = NULL, niter_max_u_grid = 5, tolerance_u_grid = 10, u_delta = 5) {
  # Using multi-block simulator with unconnected blocks
  out <- SimulateGraphical(
    n = pk, pk = n, theta = adjacency,
    implementation = SimulateAdjacency,
    topology = "random",
    nu_within = nu_within, # fully connected components by default
    nu_between = nu_between, # fully connected by default but with v_between=0
    output_matrices = output_matrices,
    v_within = v_within,
    v_between = v_between, # unconnected blocks if set to zero
    continuous = continuous,
    pd_strategy = pd_strategy,
    u = u, niter_max_u_grid = niter_max_u_grid,
    tolerance_u_grid = tolerance_u_grid, u_delta = u_delta
  )

  # Transposing the matrix for clusters of individuals
  out$data <- t(out$data)
  rownames(out$data) <- paste0("obs", 1:nrow(out$data))
  colnames(out$data) <- paste0("var", 1:ncol(out$data))

  # Updating row and column names of output matrices
  if ("omega" %in% names(out)) {
    rownames(out$omega) <- colnames(out$omega) <- rownames(out$data)
    rownames(out$phi) <- colnames(out$phi) <- rownames(out$data)
    rownames(out$C) <- colnames(out$C) <- rownames(out$data)
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

  return(out)
}


#' Simulation of predictors and associated outcome
#'
#' Simulate (i) a matrix X of n observations for \code{pk} normally distributed
#' variables, and (ii) an outcome Y obtained from a linear combination of (a
#' subset of) the \code{pk} variables in X. The outputs of this function can be
#' used to evaluate the ability of variable selection algorithms to identify,
#' among the variables in X, relevant predictors of the outcome Y.
#'
#' @inheritParams SimulateGraphical
#' @param pk number of variables in the simulated dataset X.
#' @param X matrix of predictors. With \code{X=NULL}, the matrix of predictors
#'   is simulated for \code{pk} variables and \code{n} observations. If \code{X}
#'   is provided, (a subset of) its variables will be used to simulate the
#'   outcome data.
#' @param theta optional binary vector encoding true predictors.
#' @param nu_pred density of the set of variables to be used for the simulation
#'   of the outcome data (i.e. number of true predictors over the number of
#'   potential predictors).
#' @param beta_set vector defining the (range of) coefficients used in the
#'   linear combination. If \code{continuous=FALSE}, \code{beta_set} is the set
#'   of possible coefficients. If \code{continuous=FALSE}, \code{beta_set} is
#'   the range of possible coefficients.
#' @param continuous logical indicating whether to sample coefficients from a
#'   uniform distribution between the minimum and maximum values in
#'   \code{beta_set} (\code{continuous=TRUE}) or from proposed values in
#'   \code{beta_set} (\code{continuous=FALSE}).
#' @param prop_ev proportion of explained variance in the outcome from the best
#'   linear combination of the predictors. This value must be between 0 and 1.
#' @param family type of outcome. If \code{family="gaussian"}, a normally
#'   distributed outcome is simulated. If \code{family="binomial"}, a binary
#'   outcome is simulated from a Bernouilli distribution where the probability
#'   is defined as a linear combination of (a subset of) the \code{pk}
#'   variables.
#'
#' @return A list with: \item{X}{simulated dataset. A subset of these variables
#'   were used in the linear combination defining the outcome Y.}
#'   \item{Y}{simulated outcome Y.} \item{proba}{true probability that the
#'   outcome is equal to 1. Only returned if \code{family="binomial"}.}
#'   \item{logit_proba}{logit of the true probability that the outcome is equal
#'   to 1. This is obtained as a linear combination of (a subset of) the
#'   variables in X. Only returned if \code{family="binomial"}.}
#'   \item{theta}{binary vector indicating which variables in X were used in the
#'   linear combination to obtain the outcome Y, i.e. indicating which variables
#'   in X are signal (if equal to 1) or noise (if equal to 0) variables in
#'   association with the outcome Y.} \item{beta}{true coefficients used in the
#'   linear model for simulation of the outcome Y.}
#'
#' @family simulation functions
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation (continuous outcome)
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 100, family = "gaussian", prop_ev = 0.8)
#' plot(simul$Y_pred, simul$Y) # true linear combination vs simulated outcome
#' simul <- SimulateRegression(n = 200, pk = 100, family = "gaussian", prop_ev = 0.3)
#' plot(simul$Y_pred, simul$Y) # larger residual error
#'
#' # Outcome simulation from user-provided predictor data
#' set.seed(1)
#' X <- matrix(rnorm(n = 200), nrow = 20, ncol = 10)
#' simul <- SimulateRegression(X = X, family = "gaussian")
#'
#' # Manually chosen set of predictors
#' mytheta <- c(rep(0, 8), rep(1, 2))
#' simul <- SimulateRegression(X = X, theta = mytheta, family = "gaussian")
#'
#' # Data simulation (binary outcome)
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 100, family = "binomial")
#' boxplot(simul$logit_proba ~ simul$Y) # true logit probability by simulated binary outcome
#' }
#' @export
SimulateRegression <- function(n = 100, pk = 10,
                               X = NULL, theta = NULL, adjacency = NULL,
                               nu_pred = 0.2, beta_set = c(-1, 1), continuous = FALSE,
                               prop_ev = 0.8, family = "gaussian") {
  # Checking the inputs
  if ((!is.null(X)) & (!is.null(theta))) {
    if (ncol(X) != length(theta)) {
      stop("Arguments 'X' and 'theta' are not consistent. The length of vector 'theta' must be equal to the number of columns in 'X'.")
    }
  }

  # Definition of the number of (noise+signal) predictor variables
  if (is.null(theta)) {
    p <- sum(pk)
  } else {
    pk <- p <- length(theta)
  }

  # Simulation of the predictors if not provided
  if (is.null(X)) {
    p <- sum(pk)
    X <- NULL
    for (k in 1:p) {
      X <- cbind(X, stats::rnorm(n, mean = 0, sd = 1))
    }
    X <- scale(X)
  } else {
    if (!is.null(X)) {
      n <- nrow(X)
      p <- ncol(X)
    }
  }

  # Setting column names for predictors
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("var", 1:ncol(X))
  }

  # Setting row names for predictors
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("obs", 1:nrow(X))
  }

  # Getting the binary vector of true predictors
  if (is.null(theta)) {
    theta <- stats::rbinom(p, size = 1, prob = nu_pred)
  }
  names(theta) <- colnames(X)

  # Simulating a vector of betas
  if (continuous) {
    beta <- stats::runif(p, min = min(beta_set), max = max(beta_set))
  } else {
    beta <- base::sample(beta_set, size = p, replace = TRUE)
  }
  beta <- beta * theta

  # Computing the predicted values of Y
  Y_pred <- X %*% beta

  # Estimating variance of error term to reach specified proportion of explained variance
  sd_pred_error <- sqrt(stats::var(Y_pred) * (1 / prop_ev - 1))

  # Introducing some centered gaussian error
  Y <- Y_pred + stats::rnorm(n, mean = 0, sd = sd_pred_error)

  # Compute binary outcome for logistic regression
  if (family == "binomial") {
    proba <- 1 / (1 + exp(-Y)) # inverse logit
    Y_bin <- base::cbind(stats::rbinom(n, size = 1, prob = proba))
  }

  # Return the simulated X and Y
  if (family == "binomial") {
    out <- list(X = X, Y = Y_bin, proba = proba, logit_proba = Y, logit_proba_pred = Y_pred, theta = theta, beta = beta)
  } else {
    out <- list(X = X, Y = Y, Y_pred = Y_pred, theta = theta, beta = beta)
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
#' @param nu expected density of the graph. This argument is only used for
#'   \code{topology="random"} or \code{topology="cluster"} (see argument
#'   \code{prob} in \code{\link[huge]{huge.generator}}).
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
#' adjacency <- SimulateAdjacency(pk = 20, topology = "scale-free")
#' plot(Graph(adjacency))
#' }
#' @export
SimulateAdjacency <- function(pk = 10, topology = "random", nu = 0.1, ...) {
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
#' @export
SimulateBlockAdjacency <- function(pk = 10,
                                   implementation = SimulateAdjacency, topology = "random",
                                   nu_within = 0.1, nu_between = NULL, ...) {
  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Simulation of the adjacency matrix
  if (length(pk) > 1) {
    # Initialising theta
    theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
    theta_vect <- theta[upper.tri(theta)]

    # Allowing for different densities in within and between blocks
    theta_w <- do.call(implementation, args = list(pk = pk, topology = topology, nu = nu_within, ...))
    theta_w_vect <- theta_w[upper.tri(theta_w)]
    theta_b <- do.call(implementation, args = list(pk = pk, topology = topology, nu = nu_between, ...))
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
    theta <- do.call(implementation, args = list(pk = pk, topology = topology, nu = nu_within, ...))
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
#' @examples
#' \dontrun{
#'
#' # Simulating a symmetric with 2 blocks
#' mat <- SimulateSymmetricMatrix(pk = c(5, 5))
#' }
#'
#' @keywords internal
SimulateSymmetricMatrix <- function(pk = 10, v_within = c(-1, 1), v_between = c(-0.1, 0.1), continuous = FALSE) {
  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]

  # Making as factor to allow for groups with 1 variable (for clustering)
  bigblocks_vect <- factor(bigblocks_vect, levels = seq(1, max(bigblocks)))
  block_ids <- unique(as.vector(bigblocks))

  # Building v matrix
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
#' matrix indicate pairwise conditional independence.
#'
#' @inheritParams SimulateGraphical
#' @param theta binary and symmetric adjacency matrix encoding the conditional
#'   independence structure.
#'
#' @return A list with: \item{omega}{true simulated precision matrix.}
#'   \item{u}{chosen value of u.} \item{u_grid}{grid of u values.}
#'   \item{contrast_path}{contrast values obtained for the values of u listed in
#'   u_grid.}
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of an adjacency matrix
#' theta <- SimulateBlockAdjacency(pk = c(5, 5))
#' print(theta)
#'
#' # Simulation of a precision matrix
#' simul <- SimulatePrecision(theta = theta)
#' print(simul$omega)
#' }
#' @export
SimulatePrecision <- function(pk = 10, theta,
                              v_within = c(-1, 1), v_between = c(-0.1, 0.1), continuous = FALSE,
                              pd_strategy = "diagonally_dominant",
                              u = NULL, niter_max_u_grid = 5, tolerance_u_grid = 10, u_delta = 5) {
  # Defining grid of u values if not provided
  if (is.null(u)) {
    u <- 10^-(seq(0, 5, by = 0.1))
    refining_u_grid <- TRUE
    niter_max <- 5
    tolerance <- 10
  } else {
    refining_u_grid <- FALSE
  }

  # Building v matrix
  v <- SimulateSymmetricMatrix(pk = pk, v_within = v_within, v_between = v_between, continuous = continuous)

  # Filling off-diagonal entries of the precision matrix
  omega <- theta * v

  # Calibrate u based on contrasts of the correlation matrix
  contrast <- NULL
  for (u_value in u) {
    omega_tmp <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
    C <- stats::cov2cor(solve(omega_tmp))
    contrast <- c(contrast, Contrast(C))
  }

  # Avoiding extreme values in u grid if not provided by the user
  if (refining_u_grid) {
    stop <- 0
    niter <- 1
    while (stop == 0) {
      niter <- niter + 1
      if (niter == niter_max_u_grid) {
        stop <- 1
      }
      # Satisfied with calibrated u if the argmax is not too close to the boundaries (as defined from tolerance_u_grid)
      if (any(which(contrast == max(contrast)) %in% seq(tolerance_u_grid, length(u) - tolerance_u_grid) == TRUE)) {
        stop <- 1
      } else {
        # Adding smaller values of u
        if (any(which(contrast == max(contrast)) %in% seq(1, tolerance_u_grid) == TRUE)) {
          u <- c(u, 10^-seq(min(-log10(u)) - u_delta, min(-log10(u)), by = 0.1))
        }

        # Adding larger values of u
        if (any(which(contrast == max(contrast)) %in% seq(length(u) - tolerance_u_grid, length(u)) == TRUE)) {
          u <- c(u, 10^-seq(max(-log10(u)), max(-log10(u) + u_delta), by = 0.1))
        }

        # Sorting values in u
        u <- sort(u, decreasing = TRUE)

        # Computing the contrast for all visited values of u
        contrast <- NULL
        for (u_value in u) {
          omega_tmp <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
          C <- stats::cov2cor(solve(omega_tmp))
          contrast <- c(contrast, Contrast(C))
        }
      }
    }
  }

  # Computing calibrated precision matrix
  if (length(u) > 1) {
    u_value <- u[length(contrast) - which.max(rev(contrast)) + 1] # adding smallest possible u value to the diagonal
    omega <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
  } else {
    omega <- omega_tmp
  }

  return(list(
    omega = omega,
    u = u_value, u_grid = u, contrast_path = contrast
  ))
}

#' Making positive definite
#'
#' Determines the diagonal entries of a symmetric matrix to ensure it is
#' positive definite. For this, diagonal entries of the matrix are (i) defined
#' to be higher than the sum of entries on the corresponding rows, which ensure
#' it is diagonally dominant, or (ii) set to the absolute value of the smallest
#' eigenvalue.
#'
#' @inheritParams SimulateGraphical
#' @param omega input matrix.
#' @param u_value numeric value for constant u.
#'
#' @return A positive definite matrix.
#'
#' @examples
#' \dontrun{
#'
#' # Simulation of a symmetric matrix
#' p <- 5
#' set.seed(1)
#' sigma <- matrix(rnorm(p * p), ncol = p)
#' sigma <- sigma + t(sigma)
#'
#' # Diagonal dominance
#' sigma_pd <- MakePositiveDefinite(sigma, pd = "diagonally_dominant")
#'
#' # Non-negative eigenvalues
#' sigma_pd <- MakePositiveDefinite(sigma, pd = "nonnegative_eigenvalues")
#' eigen(sigma_pd)$values
#' }
#' @export
MakePositiveDefinite <- function(omega, u_value = 0.1, pd_strategy = "diagonally_dominant") {
  # Adding a small number (u) to the diagonal
  if (pd_strategy == "diagonally_dominant") {
    diag(omega) <- apply(abs(omega), 1, sum) + u_value
  }

  # Ensuring positive eigenvalues
  if (pd_strategy == "nonnegative_eigenvalues") {
    diag(omega) <- abs(min(eigen(omega)$values)) + u_value # used in huge
  }

  return(omega)
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

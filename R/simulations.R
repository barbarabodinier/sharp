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
#' @param implementation function for simulation of the graph. By default,
#'   functionalities implemented in \code{\link[huge]{huge.generator}} are used.
#'   Alternatively, a user-defined function can be used. It must take \code{pk},
#'   \code{topology} and \code{nu} as arguments and return a
#'   \code{(sum(pk)*(sum(pk)))} binary and symmetric matrix for which diagonal
#'   entries are all equal to zero.
#' @param topology topology of the simulated graph. If using
#'   \code{implementation="huge"}, possible values are listed for the argument
#'   \code{graph} of \code{\link[huge]{huge.generator}}. These are: "random",
#'   "hub", "cluster", "band" and "scale-free".
#' @param nu expected density of the graph. If \code{implementation="huge"},
#'   this argument is only used for \code{topology="random"} or
#'   \code{topology="cluster"} (see argument \code{prob} in
#'   \code{\link[huge]{huge.generator}}).
#' @param output_matrices logical indicating if the true precision and (partial)
#'   correlation matrices should be included in the output.
#' @param v_within multiplicative factor used for diagonal blocks in simulation
#'   of the precision matrix.
#' @param v_between multiplicative factor used for off-diagonal blocks in
#'   simulation of the precision matrix.
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
#' @examples
#' \dontrun{
#'
#' # Simulation of random graph with 50 nodes
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 50, topology = "random", nu = 0.05)
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
#' simul <- SimulateGraphical(n = 100, pk = pk, nu = 0.05)
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
#' # Using user-defined function for graph simulation
#' CentralNode <- function(pk, topology = NULL, nu = NULL, hub = 1) {
#'   theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
#'   theta[hub, ] <- 1
#'   theta[, hub] <- 1
#'   diag(theta) <- 0
#'   return(theta)
#' }
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = "CentralNode")
#' plot(Graph(simul$theta)) # star
#' simul <- SimulateGraphical(n = 100, pk = 10, implementation = "CentralNode", hub = 2)
#' plot(Graph(simul$theta)) # variable 2 is the central node
#' }
#' @export
SimulateGraphical <- function(n = 100, pk = 10, implementation = SimulateAdjacency, topology = "random", nu = 0.1,
                              output_matrices = FALSE,
                              v_within = 1, v_between = 0.1,
                              pd_strategy = "diagonally_dominant",
                              u = NULL, niter_max_u_grid = 5, tolerance_u_grid = 10, u_delta = 5, ...) {
  # Defining grid of u values if not provided
  if (is.null(u)) {
    u <- 10^-(seq(0, 5, by = 0.1))
    refining_u_grid <- TRUE
    niter_max <- 5
    tolerance <- 10
  } else {
    refining_u_grid <- FALSE
  }

  # Defining number of nodes
  p <- sum(pk)

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  block_ids <- unique(as.vector(bigblocks))
  names(N_blocks) <- block_ids
  nblocks <- max(block_ids)

  # Building v matrix
  v_list <- rep(NA, length(block_ids))
  v_list[unique(diag(bigblocks))] <- v_within
  v_list[is.na(v_list)] <- v_between
  v <- bigblocks
  for (k in block_ids) {
    v[bigblocks == k] <- v_list[k]
  }

  # Simulation of the adjacency matrix
  theta <- do.call(implementation, args = list(pk = pk, topology = topology, nu = nu, ...))

  # Ensuring that there is no self-loops
  diag(theta) <- 0

  # Setting variable names
  colnames(theta) <- rownames(theta) <- paste0("var", 1:ncol(theta))

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
    u_value <- u[which.max(contrast)]
    omega <- MakePositiveDefinite(omega = omega, u_value = u_value, pd_strategy = pd_strategy)
  } else {
    omega <- omega_tmp
  }

  # Simulating the sign of the correlations (uniform)
  sign_mat <- matrix(0, nrow = nrow(omega), ncol = ncol(omega))
  sign_mat[upper.tri(sign_mat)] <- sample(c(-1, 1), size = sum(upper.tri(omega)), replace = TRUE)
  sign_mat <- sign_mat + t(sign_mat)
  diag(sign_mat) <- 1
  omega <- omega * sign_mat

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
      u = u_value, u_grid = u, contrast_path = contrast
    ))
  } else {
    return(list(data = x, theta = theta))
  }
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
#' @param sd_pred_error residual standard deviation. Used only if
#'   \code{family="gaussian"}.
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
#'   \item{theta_pred}{binary vector indicating which variables in X were used
#'   in the linear combination to obtain the outcome Y, i.e. indicating which
#'   variables in X are signal (if equal to 1) or noise (if equal to 0)
#'   variables in association with the outcome Y.} \item{beta}{true coefficients
#'   used in the linear model for simulation of the outcome Y.}
#'
#' @family simulation functions
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation (continuous outcome)
#' simul <- SimulateRegression(n = 200, pk = 100, family = "gaussian")
#' plot(simul$Y_pred, simul$Y) # true linear combination vs simulated outcome
#' simul <- SimulateRegression(n = 200, pk = 100, family = "gaussian", sd_pred_error = 5)
#' plot(simul$Y_pred, simul$Y) # larger residual error
#'
#' # Data simulation (binary outcome)
#' simul <- SimulateRegression(n = 200, pk = 100, family = "binomial")
#' boxplot(simul$logit_proba ~ simul$Y) # true logit probability by simulated binary outcome
#' }
#' @export
SimulateRegression <- function(n = 100, pk = 10, X = NULL, nu_pred = 0.2,
                               beta_set = c(-1, 1), continuous = FALSE,
                               sd_pred_error = 1, family = "gaussian") {
  # Simulation of the predictors
  if (is.null(X)) {
    p <- sum(pk)
    X <- NULL
    for (k in 1:p) {
      X <- cbind(X, stats::rnorm(n, mean = 0, sd = 1))
    }
    X <- scale(X)
  } else {
    n <- nrow(X)
    p <- ncol(X)
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
  theta_pred <- stats::rbinom(p, size = 1, prob = nu_pred)
  names(theta_pred) <- colnames(X)

  # Simulating a vector of betas
  if (continuous) {
    beta <- stats::runif(p, min = min(beta_set), max = max(beta_set))
  } else {
    beta <- base::sample(beta_set, size = p, replace = TRUE)
  }
  beta <- beta * theta_pred

  # Computing the predicted values of Y
  Y_pred <- X %*% beta

  # Introducing some centered gaussian error
  Y <- Y_pred + stats::rnorm(n, mean = 0, sd = sd_pred_error)

  # Compute binary outcome for logistic regression
  if (family == "binomial") {
    proba <- 1 / (1 + exp(-Y)) # inverse logit
    Y_bin <- base::cbind(stats::rbinom(n, size = 1, prob = proba))
  }

  # Return the simulated X and Y
  if (family == "binomial") {
    out <- list(X = X, Y = Y_bin, proba = proba, logit_proba = Y, logit_proba_pred = Y_pred, theta_pred = theta_pred, beta = beta)
  } else {
    out <- list(X = X, Y = Y, Y_pred = Y_pred, theta_pred = theta_pred, beta = beta)
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
#' Computes the matrix contrast, defined as the number of visited bins of
#' entries with a specified number of digits.
#'
#' @param mat input matrix.
#' @param digits number of digits to use.
#'
#' @return A single number, the contrast of the input matrix.
Contrast <- function(mat, digits = 3) {
  return(length(unique(round(as.vector(abs(mat)), digits = digits))))
}

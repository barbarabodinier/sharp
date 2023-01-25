#' Merging stability selection outputs
#'
#' Merges the outputs from two runs of \code{\link{VariableSelection}},
#' \code{\link{GraphicalModel}} or \code{\link{Clustering}}. The two runs must
#' have been done using the same \code{methods} and the same \code{params} but
#' with different \code{seed}s. The combined output will contain results based
#' on iterations from both \code{stability1} and \code{stability2}. This
#' function can be used for parallelisation.
#'
#' @param stability1 output from a first run of \code{\link{VariableSelection}},
#'   \code{\link{GraphicalModel}}, or \code{\link{Clustering}}.
#' @param stability2 output from a second run of
#'   \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}, or
#'   \code{\link{Clustering}}.
#' @param include_beta logical indicating if the beta coefficients of visited
#'   models should be concatenated. Only applicable to variable selection or
#'   clustering.
#'
#' @return A single output of the same format.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' \donttest{
#' ## Variable selection
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'
#' # Two runs
#' stab1 <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, seed = 1, K = 10)
#' stab2 <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2, include_beta = FALSE)
#' str(stab)
#'
#'
#' ## Graphical modelling
#'
#' # Data simulation
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Two runs
#' stab1 <- GraphicalModel(xdata = simul$data, seed = 1, K = 10)
#' stab2 <- GraphicalModel(xdata = simul$data, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2)
#' str(stab)
#'
#'
#' ## Clustering
#'
#' # Data simulation
#' simul <- SimulateClustering(n = c(15, 15, 15))
#'
#' # Two runs
#' stab1 <- Clustering(xdata = simul$data, seed = 1)
#' stab2 <- Clustering(xdata = simul$data, seed = 2)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2)
#' str(stab)
#' }
#' @export
Combine <- function(stability1, stability2, include_beta = TRUE) {
  if (!inherits(stability1, c("graphical_model", "variable_selection", "clustering"))) {
    stop("Invalid inputs. This function only applies to outputs from GraphicalModel(), VariableSelection() or Clustering().")
  }
  if (class(stability1) != class(stability2)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They must be generated from the same function.")
  }
  if (!all(is.na(stability1$Lambda))) {
    if (any(stability1$Lambda != stability2$Lambda)) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different Lambdas.")
    }
  }
  if (any(do.call(c, stability1$methods) != do.call(c, stability2$methods))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (inherits(stability1, "clustering")) {
    if (any(do.call(c, stability1$params[c("pk", "n", "tau")]) != do.call(c, stability2$params[c("pk", "n", "tau")]))) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
    }
  } else {
    if (any(do.call(c, stability1$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]) != do.call(c, stability2$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]))) {
      stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
    }
  }
  if (stability1$params$seed == stability2$params$seed) {
    stop("Arguments 'stability1' and 'stability2' were obtained using the same seed.")
  }
  if (any(stability1$sign != stability2$sign)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }

  # Identifying the type of model (variable/pairs)
  if (inherits(stability1, "graphical_model")) {
    graph <- TRUE
  } else {
    graph <- FALSE
  }

  # Extracting the parameters (some are NULL depending on the method)
  nc <- stability1$nc
  Lambda <- stability1$Lambda
  mysign <- stability1$sign
  pk <- stability1$params$pk
  pi_list <- stability1$params$pi_list
  n_cat <- stability1$params$n_cat
  Sequential_template <- stability1$params$Sequential_template
  PFER_method <- stability1$methods$PFER_method
  PFER_thr <- stability1$params$PFER_thr
  FDP_thr <- stability1$params$FDP_thr
  mymethods <- stability1$methods
  myparams <- stability1$params

  # Creating matrix with block indices (specific to graphical models)
  nblocks <- 1
  N_blocks <- stability1$params$pk
  names(N_blocks) <- 1
  if (stability1$methods$type == "graphical_model") { # to avoid memory issues in high dimensional variable selection
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- unique(as.vector(bigblocks_vect))
    names(N_blocks) <- blocks
    nblocks <- max(blocks)
  }

  # Preparing the PFER and FDP thresholds
  if (length(PFER_thr) == 1) {
    PFER_thr_blocks <- ceiling(prop.table(N_blocks) * PFER_thr)
  } else {
    if (length(PFER_thr) == nblocks) {
      PFER_thr_blocks <- PFER_thr
    }
  }
  if (length(FDP_thr) == 1) {
    FDP_thr_blocks <- rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr) == nblocks) {
      FDP_thr_blocks <- FDP_thr
    }
  }

  # Computing the total number of iterations
  K <- stability1$params$K + stability2$params$K
  myparams$K <- K
  myparams$seed <- Inf

  # Computing selection proportions
  bigstab <- NULL
  if (!is.null(stability1$selprop)) {
    bigstab <- array(NA, dim = dim(stability1$selprop), dimnames = dimnames(stability1$selprop))
    if (graph) {
      for (k in 1:dim(stability1$selprop)[3]) {
        bigstab[, , k] <- (stability1$selprop[, , k] * stability1$params$K + stability2$selprop[, , k] * stability2$params$K) / K
      }
    } else {
      for (k in 1:nrow(stability1$selprop)) {
        bigstab[k, ] <- (stability1$selprop[k, ] * stability1$params$K + stability2$selprop[k, ] * stability2$params$K) / K
      }
    }
  }

  # Computing co-membership proportions
  if (inherits(stability1, "clustering")) {
    coprop <- array(NA, dim = dim(stability1$coprop), dimnames = dimnames(stability1$coprop))
    for (k in 1:dim(stability1$coprop)[3]) {
      coprop[, , k] <- (stability1$coprop[, , k] * stability1$params$K + stability2$coprop[, , k] * stability2$params$K) / K
    }
  }

  # Concatenating the beta coefficients
  if (inherits(stability1, c("variable_selection", "clustering"))) {
    if (include_beta) {
      Beta <- NULL
      if (!is.null(stability1$Beta)) {
        if (length(dim(stability1$Beta)) == 4) {
          Beta <- array(NA,
            dim = c(dim(stability1$Beta)[1:2], dim(stability1$Beta)[3] + dim(stability2$Beta)[3], dim(stability1$Beta)[4]),
            dimnames = list(
              dimnames(stability1$Beta)[[1]], dimnames(stability1$Beta)[[2]],
              c(dimnames(stability1$Beta)[[3]], dimnames(stability2$Beta)[[3]]), dimnames(stability1$Beta)[[4]]
            )
          )
          Beta[, , 1:dim(stability1$Beta)[3], ] <- stability1$Beta
          Beta[, , (dim(stability1$Beta)[3] + 1):dim(Beta)[3], ] <- stability2$Beta
        } else {
          Beta <- array(NA,
            dim = c(dim(stability1$Beta)[1:2], dim(stability1$Beta)[3] + dim(stability2$Beta)[3]),
            dimnames = list(
              dimnames(stability1$Beta)[[1]], dimnames(stability1$Beta)[[2]],
              c(dimnames(stability1$Beta)[[3]], dimnames(stability2$Beta)[[3]])
            )
          )
          Beta[, , 1:dim(stability1$Beta)[3]] <- stability1$Beta
          Beta[, , (dim(stability1$Beta)[3] + 1):dim(Beta)[3]] <- stability2$Beta
        }
      }
    }
  }

  # Computation of the stability score
  if (inherits(stability1, "graphical_model")) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  }
  if (inherits(stability1, "variable_selection")) {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  }
  if (inherits(stability1, "clustering")) {
    Sc <- matrix(NA, nrow = dim(coprop)[3], ncol = 1)
    for (k in 1:dim(coprop)[3]) {
      Sc[k, 1] <- ConsensusScore(coprop = coprop[, , k], nc = nc[k], K = K, linkage = stability1$methods$linkage)
    }
    Q <- 1 / K * (stability1$params$K * stability1$Q + stability2$params$K * stability2$Q) # weighted average
  }

  # Preparing output
  if (inherits(stability1, "graphical_model")) {
    if (nblocks == 1) {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      )
    }
  }

  if (inherits(stability1, "variable_selection")) {
    if (include_beta) {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        Beta = Beta,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        methods = mymethods,
        params = myparams
      )
    }
  }

  if (inherits(stability1, "clustering")) {
    if (include_beta) {
      out <- list(
        Sc = Sc,
        nc = nc,
        Lambda = Lambda,
        Q = Q,
        coprop = coprop,
        Beta = Beta,
        selprop = bigstab,
        methods = mymethods,
        params = myparams
      )
    } else {
      out <- list(
        Sc = Sc,
        nc = nc,
        Lambda = Lambda,
        Q = Q,
        coprop = coprop,
        selprop = bigstab,
        methods = mymethods,
        params = myparams
      )
    }
  }

  # Defining the class
  class(out) <- class(stability1)

  return(out)
}

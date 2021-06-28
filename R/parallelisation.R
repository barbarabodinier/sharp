#' Merging stability selection outputs
#'
#' Merges the outputs from \code{\link{VariableSelection}} or
#' \code{\link{GraphicalModel}}. This function can be used for parallelisation.
#'
#' @param stability1 output from a first run of \code{\link{VariableSelection}}
#'   or \code{\link{GraphicalModel}}.
#' @param stability2 output from a second run of \code{\link{VariableSelection}}
#'   or \code{\link{GraphicalModel}}.
#' @param include_beta logical indicating if the beta coefficients of visited
#'   models should be concatenated.
#'
#' @return A single output of the same format.
#'
#' @details The two runs must have been done using the same data, the same grids
#'   of parameters, the same methods (arguments \code{implementation},
#'   \code{resampling}, \code{PFER_method}, as well as \code{start} for
#'   graphical models), the same \code{tau} and \code{pk} (graphical models
#'   only), and the same thresholds in PFER and FDP, but with different seeds.
#'   The combined output will be the equivalent of running the model for which
#'   the number of iterations is the sum of the number of iterations used for
#'   \code{stability1} and \code{stability2}.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#'
#' ## Variable selection
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'
#' # Two runs
#' stab1 <- VariableSelection(xdata = simul$X, ydata = simul$Y, seed = 1, K = 10)
#' stab2 <- VariableSelection(xdata = simul$X, ydata = simul$Y, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2, include_beta = FALSE)
#' print(stab$params$K)
#'
#' ## Graphical modelling
#' # Data simulation
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Two runs
#' stab1 <- GraphicalModel(xdata = simul$data, seed = 1, K = 10)
#' stab2 <- GraphicalModel(xdata = simul$data, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2)
#' print(stab$params$K)
#'
#' ## Clustering
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(
#'   n = c(5, 5, 5), pk = 100,
#'   v_within = c(-1, -0.5), continuous = TRUE
#' )
#'
#' # Two runs
#' stab1 <- Clustering(xdata = simul$data, seed = 1, K = 10)
#' stab2 <- Clustering(xdata = simul$data, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2)
#' print(stab$params$K)
#' @export
Combine <- function(stability1, stability2, include_beta = TRUE) {
  if (any(stability1$Lambda != stability2$Lambda)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different Lambdas.")
  }
  if (any(do.call(c, stability1$methods) != do.call(c, stability2$methods))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (any(do.call(c, stability1$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]) != do.call(c, stability2$params[c("pk", "n", "tau", "PFER_thr", "FDP_thr")]))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (stability1$params$seed == stability2$params$seed) {
    warning("Arguments 'stability1' and 'stability2' were obtained using the same seed.")
  }
  if (any(stability1$sign != stability2$sign)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }

  # Identifying the type of model (variable/pairs)
  if (stability1$methods$type %in% c("graphical_model", "clustering")) {
    graph <- TRUE
  } else {
    graph <- FALSE
  }

  # Extracting the parameters
  Lambda <- stability1$Lambda
  if (stability1$methods$type == "graphical_model") {
    mysign <- stability1$sign
  }
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

  # Computing selection propotions
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

  # Concatenating the beta coefficients
  if (stability1$methods$type == "variable_selection") {
    if (include_beta) {
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

  # Computation of the stability score
  if (stability1$methods$type == "graphical_model") {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  } else {
    metrics <- StabilityMetrics(
      selprop = bigstab, pk = NULL, pi_list = pi_list, K = K, n_cat = n_cat,
      Sequential_template = Sequential_template, graph = graph,
      PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
    )
  }

  # Preparing output
  if (stability1$methods$type == "graphical_model") {
    if (nblocks == 1) {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      ))
    } else {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d,
        selprop = bigstab,
        sign = mysign,
        methods = mymethods,
        params = myparams
      ))
    }
  }

  if (stability1$methods$type == "variable_selection") {
    if (include_beta) {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        Beta = Beta,
        methods = mymethods,
        params = myparams
      ))
    } else {
      return(list(
        S = metrics$S, Lambda = Lambda,
        Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
        PFER = metrics$PFER, FDP = metrics$FDP,
        S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
        selprop = bigstab,
        methods = mymethods,
        params = myparams
      ))
    }
  }

  if (stability1$methods$type == "clustering") {
    return(list(
      S = metrics$S, Lambda = Lambda,
      Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
      PFER = metrics$PFER, FDP = metrics$FDP,
      S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
      selprop = bigstab,
      methods = mymethods,
      params = myparams
    ))
  }
}

#' Merging stability selection outputs
#'
#' Merges the outputs from \code{\link{VariableSelection}} or
#' \code{\link{GraphicalModel}}. This can be useful for parallelisation.
#'
#' @param stability1 output from a first run of \code{\link{VariableSelection}}
#'   or \code{\link{GraphicalModel}}.
#' @param stability2 output from a second run of \code{\link{VariableSelection}}
#'   or \code{\link{GraphicalModel}}.
#' @param graph logical indicating if \code{\link{GraphicalModel}} (graph=TRUE)
#'   or \code{\link{VariableSelection}} (graph=FALSE) were used to generate the
#'   "stability1" and "stability2".
#'
#' @return a single output with the same format.
#'
#' @details The two runs must have been done using the same data, the same grids
#' of parameters, the same methods (arguments "implementation", "resampling",
#' "PFER_method", as well as "start" for graphical models), the same "tau" and
#' "pk" (graphical models only), and the same thresholds in PFER and FDP, but
#' with different seeds. The combined output will be the equivalent of running
#' the model for which the number of iterations is the sum of the number of
#' iterations used for "stability1" and "stability2".
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{GraphicalModel}}
#'
#' @examples
#' # Data simulation
#' simul <- SimulateGraphical(pk = 20)
#'
#' # Two runs
#' stab1 <- GraphicalModel(data = simul$data, seed = 1, K = 10)
#' stab2 <- GraphicalModel(data = simul$data, seed = 2, K = 10)
#'
#' # Merging the outputs
#' stab <- Combine(stability1 = stab1, stability2 = stab2, graph = TRUE)
#' print(stab$params$K)
#' @export
Combine <- function(stability1, stability2, graph = TRUE) {
  if (any(stability1$Lambda != stability2$Lambda)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different Lambdas.")
  }
  if (any(do.call(c, stability1$methods) != do.call(c, stability2$methods))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (any(do.call(c, stability1$params[c("pk", "tau", "PFER_thr", "FDP_thr")]) != do.call(c, stability2$params[c("pk", "tau", "PFER_thr", "FDP_thr")]))) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (stability1$params$seed == stability2$params$seed) {
    warning("Arguments 'stability1' and 'stability2' were obtained using the same seed.")
  }
  if (any(stability1$params$data != stability2$params$data)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$params$xdata != stability2$params$xdata)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$params$ydata != stability2$params$ydata)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$sign != stability2$sign)) {
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }

  # Extracting the parameters
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

  # Creating matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

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

  # Computation of the stability score
  metrics <- StabilityMetrics(
    bigstab = bigstab, pk = pk, pi_list = pi_list, K = K, n_cat = n_cat,
    Sequential_template = Sequential_template, graph = graph,
    PFER_method = PFER_method, PFER_thr_blocks = PFER_thr_blocks, FDP_thr_blocks = FDP_thr_blocks
  )

  if (nblocks == 1) {
    return(list(
      S = metrics$S, Lambda = Lambda,
      Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
      PFER = metrics$PFER, FDP = metrics$FDP,
      S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, FDP_2d = metrics$FDP_2d,
      selprop = bigstab, sign = mysign,
      methods = mymethods,
      params = myparams
    ))
  } else {
    return(list(
      S = metrics$S, Lambda = Lambda,
      Q = metrics$Q, Q_s = metrics$Q_s, P = metrics$P,
      PFER = metrics$PFER, FDP = metrics$FDP,
      S_2d = metrics$S_2d,
      selprop = bigstab, sign = mysign,
      methods = mymethods,
      params = myparams
    ))
  }
}

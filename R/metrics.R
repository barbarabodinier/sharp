#' Stability selection metrics
#'
#' Computes the stability score (see \code{\link{StabilityScore}}) and
#' upper-bounds of the \code{\link{PFER}} and \code{\link{FDP}} from selection
#' proportions of models with a given parameter controlling the sparsity of the
#' underlying algorithm and for different thresholds in selection proportions.
#'
#' @inheritParams GraphicalModel
#' @param selprop array of selection proportions.
#' @param PFER_thr_blocks vector of block-specific thresholds in PFER for
#'   constrained calibration by error control. If \code{PFER_thr=Inf} and
#'   \code{FDP_thr=Inf}, unconstrained calibration is used.
#' @param FDP_thr_blocks vector of block-specific thresholds in the expected
#'   proportion of falsely selected features (or False Discovery Proportion,
#'   FDP) for constrained calibration by error control. If \code{PFER_thr=Inf}
#'   and \code{FDP_thr=Inf}, unconstrained calibration is used.
#' @param Sequential_template logical matrix encoding the type of procedure to
#'   use for data with multiple blocks in stability selection graphical
#'   modelling. For multi-block estimation, the stability selection model is
#'   constructed as the union of block-specific stable edges estimated while the
#'   others are weakly penalised (\code{TRUE} only for the block currently being
#'   calibrated and \code{FALSE} for other blocks). Other approaches with joint
#'   calibration of the blocks are allowed (all entries are set to \code{TRUE}).
#' @param graph logical indicating if stability selection is performed in a
#'   regression (if \code{FALSE}) or graphical (if \code{TRUE})
#'   framework.
#' @param group vector encoding the grouping structure among predictors. This
#'   argument indicates the number of variables in each group and only needs to
#'   be provided for group (but not sparse group) penalisation.
#'
#' @return A list with: \item{S}{a matrix of the best (block-specific) stability
#'   scores for different (sets of) penalty parameters. In multi-block stability
#'   selection, rows correspond to different sets of penalty parameters, (values
#'   are stored in the output "Lambda") and columns correspond to different
#'   blocks.} \item{Lambda}{a matrix of (block-specific) penalty parameters. In
#'   multi-block stability selection, rows correspond to sets of penalty
#'   parameters and columns correspond to different blocks.} \item{Q}{a matrix
#'   of average numbers of (block-specific) edges selected by the underlying
#'   algorithm for different (sets of) penalty parameters. In multi-block
#'   stability selection, rows correspond to different sets of penalty
#'   parameters, (values are stored in the output "Lambda") and columns
#'   correspond to different blocks.} \item{Q_s}{a matrix of calibrated numbers
#'   of (block-specific) stable edges for different (sets of) penalty
#'   parameters. In multi-block stability selection, rows correspond to
#'   different sets of penalty parameters, (values are stored in the output
#'   "Lambda") and columns correspond to different blocks.} \item{P}{a matrix of
#'   calibrated (block-specific) thresholds in selection proportions for
#'   different (sets of) penalty parameters. In multi-block stability selection,
#'   rows correspond to different sets of penalty parameters, (values are stored
#'   in the output "Lambda") and columns correspond to different blocks.}
#'   \item{PFER}{a matrix of computed (block-specific) upper-bounds in PFER of
#'   calibrated graphs for different (sets of) penalty parameters. In
#'   multi-block stability selection, rows correspond to different sets of
#'   penalty parameters, (values are stored in the output "Lambda") and columns
#'   correspond to different blocks.} \item{FDP}{a matrix of computed
#'   (block-specific) upper-bounds in FDP of calibrated stability selection
#'   models for different (sets of) penalty parameters. In multi-block stability
#'   selection, rows correspond to different sets of penalty parameters, (values
#'   are stored in the output "Lambda") and columns correspond to different
#'   blocks.} \item{S_2d}{an array of (block-specific) stability scores obtained
#'   with different combinations of parameters. Rows correspond to different
#'   (sets of) penalty parameters and columns correspond to different thresholds
#'   in selection proportions. In multi-block stability selection, indices along
#'   the third dimension correspond to different blocks.} \item{PFER_2d}{an
#'   array of computed upper-bounds of PFER obtained with different combinations
#'   of parameters. Rows correspond to different penalty parameters and columns
#'   correspond to different thresholds in selection proportions. Not available
#'   in multi-block stability selection graphical modelling.} \item{FDP_2d}{an
#'   array of computed upper-bounds of FDP obtained with different combinations
#'   of parameters. Rows correspond to different penalty parameters and columns
#'   correspond to different thresholds in selection proportions. Not available
#'   in multi-block stability selection graphical modelling.}
#'
#' @family stability metric functions
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#'   \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{stabilityselectionSS}{sharp}
#'
#' @examples
#' ## Sparse or sparse group penalisation
#'
#' # Simulating set of selection proportions
#' set.seed(1)
#' selprop <- matrix(round(runif(n = 20), digits = 2), nrow = 2)
#'
#' # Computing stability scores for different thresholds
#' metrics <- StabilityMetrics(
#'   selprop = selprop, pi = c(0.6, 0.7, 0.8),
#'   K = 100, graph = FALSE
#' )
#'
#'
#' ## Group penalisation
#'
#' # Simulating set of selection proportions
#' set.seed(1)
#' selprop <- matrix(round(runif(n = 6), digits = 2), nrow = 2)
#' selprop <- cbind(
#'   selprop[, 1], selprop[, 1],
#'   selprop[, 2], selprop[, 2],
#'   matrix(rep(selprop[, 3], each = 6), nrow = 2, byrow = TRUE)
#' )
#'
#' # Computing stability scores for different thresholds
#' metrics <- StabilityMetrics(
#'   selprop = selprop, pi = c(0.6, 0.7, 0.8),
#'   K = 100, graph = FALSE, group = c(2, 2, 6)
#' )
#' @export
StabilityMetrics <- function(selprop, pk = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, n_cat = 3,
                             PFER_method = "MB", PFER_thr_blocks = Inf, FDP_thr_blocks = Inf,
                             Sequential_template = NULL, graph = TRUE, group = NULL) {
  if (graph) {
    nlambda <- dim(selprop)[3]
  } else {
    nlambda <- nrow(selprop)
  }

  # Extracting pk
  if (is.null(pk)) {
    pk <- ncol(selprop)
  }

  if (is.null(Sequential_template)) {
    Sequential_template <- matrix(TRUE, nrow = nlambda, ncol = 1)
  }

  # Create matrix with block indices
  nblocks <- 1
  if (graph) { # to avoid memory issues in high dimensional variable selection
    bigblocks <- BlockMatrix(pk)
    bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
    N_blocks <- unname(table(bigblocks_vect))
    blocks <- unique(as.vector(bigblocks_vect))
    names(N_blocks) <- blocks
    nblocks <- max(blocks)
  }

  # Initialising objects to be filled
  Q <- Q_s <- P <- matrix(NA, nrow = nlambda, ncol = nblocks)
  best_loglik <- best_PFER <- best_FDP <- matrix(NA, nrow = nlambda, ncol = nblocks)
  if (nblocks == 1) {
    loglik <- PFER <- FDP <- matrix(NA, ncol = length(pi_list), nrow = nlambda)
  } else {
    loglik <- array(NA, dim = c(nlambda, length(pi_list), nblocks))
  }

  # Computing the metrics for each value of lambda
  for (k in 1:nlambda) {
    # Extracting corresponding selection proportions
    if (graph) {
      stab_iter <- selprop[, , k]
    } else {
      stab_iter <- selprop[k, ]
    }

    # Computing stability score with block-specific pi
    for (block_id in 1:nblocks) {
      if (Sequential_template[k, block_id]) {
        if (graph) {
          stab_iter_block <- stab_iter[(bigblocks == block_id) & (upper.tri(bigblocks))] # selection proportions in the block
        } else {
          stab_iter_block <- stab_iter
        }

        # Using group penalisation (extracting one per group)
        if (!is.null(group)) {
          stab_iter_block <- stab_iter_block[cumsum(group)]
        }

        q_block <- round(sum(stab_iter_block, na.rm = TRUE)) # average number of edges selected by the original procedure in the block
        Q[k, block_id] <- q_block
        N_block <- length(stab_iter_block) # maximum number of edges in the block
        tmp_loglik <- tmp_PFERs <- tmp_FDPs <- rep(NA, length(pi_list))

        # Computing error rates and stability score for different values of pi
        for (j in 1:length(pi_list)) {
          pi <- pi_list[j]
          tmp_PFERs[j] <- PFER(q = q_block, pi = pi, N = N_block, K = K, PFER_method = PFER_method)
          tmp_FDPs[j] <- FDP(selprop = stab_iter_block, PFER = tmp_PFERs[j], pi = pi)
          if ((tmp_PFERs[j] <= PFER_thr_blocks[block_id]) & (tmp_FDPs[j] <= FDP_thr_blocks[block_id])) {
            # Computing stability score (group penalisation is accounted for above so no need here)
            tmp_loglik[j] <- StabilityScore(selprop = stab_iter_block, pi_list = pi, K = K, n_cat = n_cat, group = NULL)
          }
        }

        # Storing stability score in a matrix if only one block
        if (nblocks == 1) {
          loglik[k, ] <- tmp_loglik
          PFER[k, ] <- tmp_PFERs
          FDP[k, ] <- tmp_FDPs
        } else {
          loglik[k, , block_id] <- tmp_loglik
        }

        # Keeping best stability score and other parameters at the max
        if (any(!is.na(tmp_loglik))) {
          tmp_loglik[is.na(tmp_loglik)] <- 0
          myid <- which.max(tmp_loglik)
          tmp_loglik[which(tmp_loglik == 0)] <- NA
          best_loglik[k, block_id] <- tmp_loglik[myid]
          P[k, block_id] <- pi_list[myid]
          Q_s[k, block_id] <- sum(stab_iter_block >= pi_list[myid], na.rm = TRUE)
          best_PFER[k, block_id] <- tmp_PFERs[myid]
          best_FDP[k, block_id] <- tmp_FDPs[myid]
        }
      }
    }
  }
  best_loglik_blocks <- best_loglik
  best_loglik <- matrix(apply(best_loglik, 1, sum), ncol = 1)

  if (nblocks == 1) {
    return(list(
      S = best_loglik_blocks,
      Q = Q, Q_s = Q_s, P = P,
      PFER = best_PFER, FDP = best_FDP,
      S_2d = loglik, PFER_2d = PFER, FDP_2d = FDP
    ))
  } else {
    return(list(
      S = best_loglik_blocks,
      Q = Q, Q_s = Q_s, P = P,
      PFER = best_PFER, FDP = best_FDP,
      S_2d = loglik
    ))
  }
}

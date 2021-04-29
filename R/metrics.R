#' Stability selection metrics
#'
#' This function can be used to compute the stability score and upper-bounds of
#' the PFER and FDP for stability selection models.
#'
#' @param bigstab array of selection proportions.
#' @param pk vector encoding the grouping structure. Only used for multi-block
#'   stability selection graphical models, see \code{\link{GraphicalModel}}.
#' @param pi_list grid of values for the threshold in selection proportion. With
#'   n_cat=3, these values must be between 0.5 and 1. With n_cat=2, these values
#'   must be between 0 and 1.
#' @param K number of resampling iterations.
#' @param n_cat number of categories used to compute the stability score.
#'   Possible values are 2 or 3.
#' @param PFER_method method used to compute the expected number of False
#'   Positives, (or Per Family Error Rate, PFER). With PFER_method="MB", the
#'   method proposed by Meinshausen and Bühlmann (2010) is used. With
#'   PFER_method="SS", the method proposed by Shah and Samworth (2013) under the
#'   assumption of unimodality is used.
#' @param PFER_thr_blocks (block-specific) thresholds in PFER for constrained
#'   calibration by error control. With PFER_thr_blocks=Inf and
#'   FDP_thr_blocks=Inf, unconstrained calibration is used.
#' @param FDP_thr_blocks (block-specific) thresholds in the expected proportion
#'   of falsely selected edges (or False Discovery Proportion, FDP) for
#'   constrained calibration by error control. With PFER_thr_blocks=Inf and
#'   FDP_thr_blocks=Inf, unconstrained calibration is used.
#' @param Sequential_template logical matrix encoding the type of procedure to
#'   use for data with multiple blocks in stability selection graphical
#'   modelling. For multi-block estimation, the stability selection model is
#'   constructed as the union of block-specific stable edges estimated while the
#'   others are weakly penalised (TRUE only for the block currently being
#'   calibrated and FALSE for other blocks). Other approaches with joint
#'   calibration of the blocks are allowed (all entries are set to TRUE).
#' @param graph logical indicating if stability selection is performed in a
#'   regression (FALSE) or graphical (TRUE) framework.
#'
#' @return A list with: \item{S}{a matrix of the best (block-specific) stability
#'   scores for different (sets of) penalty parameters. In multi-block stability
#'   selection, rows correspond to different sets of penalty parameters, (values
#'   are stored in the output "Lambda") and columns correspond to different
#'   blocks.} \item{Lambda}{a matrix of (block-specific) penalty parameters. In
#'   multi-block stability selection, rows correspond to sets of penalty
#'   parameters and columns correspond to different blocks.} \item{Q}{a matrix
#'   of average numbers of (block-specific) edges selected by the underlying
#'   algorihm for different (sets of) penalty parameters. In multi-block
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
#'   (sets of) penalty parameters and columns correspond to different tresholds
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
StabilityMetrics <- function(bigstab, pk = NULL, pi_list = seq(0.6, 0.9, by = 0.01),
                             K = 100, n_cat = 3,
                             PFER_method = "MB", PFER_thr_blocks = Inf, FDP_thr_blocks = Inf,
                             Sequential_template = NULL, graph = TRUE) {
  if (graph) {
    nlambda <- dim(bigstab)[3]
  } else {
    nlambda <- nrow(bigstab)
  }

  if (is.null(pk)) {
    pk <- ncol(bigstab)
  }

  if (is.null(Sequential_template)) {
    Sequential_template <- matrix(TRUE, nrow = nlambda, ncol = 1)
  }

  # Create matrix with block indices
  bigblocks <- BlockMatrix(pk)
  bigblocks_vect <- bigblocks[upper.tri(bigblocks)]
  N_blocks <- unname(table(bigblocks_vect))
  blocks <- unique(as.vector(bigblocks_vect))
  names(N_blocks) <- blocks
  nblocks <- max(blocks)

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
      stab_iter <- bigstab[, , k]
    } else {
      stab_iter <- bigstab[k, ]
    }

    # Computing stability score with block-specific pi
    for (block_id in 1:nblocks) {
      if (Sequential_template[k, block_id]) {
        if (graph) {
          stab_iter_block <- stab_iter[(bigblocks == block_id) & (upper.tri(bigblocks))] # selection proportions in the block
        } else {
          stab_iter_block <- stab_iter
        }
        q_block <- round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
        Q[k, block_id] <- q_block
        N_block <- length(stab_iter_block) # maximum number of edges in the block
        tmp_loglik <- tmp_PFERs <- tmp_FDPs <- rep(NA, length(pi_list))

        # Computing error rates and stability score for different values of pi
        for (j in 1:length(pi_list)) {
          pi <- pi_list[j]
          tmp_PFERs[j] <- PFER(q = q_block, pi = pi, N = N_block, K = K, PFER_method = PFER_method)
          tmp_FDPs[j] <- FDP(PFER = tmp_PFERs[j], pi = pi, stab_iter = stab_iter_block)
          if ((tmp_PFERs[j] <= PFER_thr_blocks[block_id]) & (tmp_FDPs[j] <= FDP_thr_blocks[block_id])) {
            tmp_loglik[j] <- StabilityScore(stab_iter = stab_iter_block, pi = pi, K = K, n_cat = n_cat)
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
          Q_s[k, block_id] <- sum(stab_iter_block >= pi_list[myid])
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

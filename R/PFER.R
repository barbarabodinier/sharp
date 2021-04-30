#' Per Family Error Rate
#'
#' Computes the upper-bound of the PFER of a stability selection model using the
#' methods proposed by Meinshausen and Bühlmann (2010) or Shah and Samworth
#' (2013). In stability selection, the PFER corresponds to the expected number
#' of stably selected features that are not relevant to the outcome (i.e. False
#' Positives).
#'
#' @param q average number of features selected by the underlying algorithm.
#' @param pi vector of thresholds in selection proportions.
#' @param N total number of features.
#' @param K number of resampling iterations.
#' @param PFER_method method used to compute the expected number of False
#'   Positives, (or Per Family Error Rate, PFER). With PFER_method="MB", the
#'   method proposed by Meinshausen and Bühlmann (2010) is used. With
#'   PFER_method="SS", the method proposed by Shah and Samworth (2013) under the
#'   assumption of unimodality is used.
#'
#' @return The upper-bound of the PFER.
#'
#' @references \insertRef{stabilityselectionMB}{focus}
#'
#' \insertRef{stabilityselectionSS}{focus}
#'
#' @family stability metric functions
#'
#' @examples
#'
#' # Example: 10 out of 50 features selected on average by underlying algorithm
#' # for stability selection with 100 iterations and threshold of 0.8
#' pfer_mb <- PFER(q = 10, pi = 0.8, N = 50, K = 100, PFER_method = "MB")
#' pfer_ss <- PFER(q = 10, pi = 0.8, N = 50, K = 100, PFER_method = "SS")
#' @export
PFER <- function(q, pi, N, K, PFER_method = "MB") {
  # Checking the inputs (PFER_method)
  PFER_method <- as.character(PFER_method)
  if ((length(PFER_method) != 1) | (!PFER_method %in% c("MB", "SS"))) {
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  if (pi > 0.5) {
    # Computing upper-bound of the PFER using approach proposed by MB
    if (PFER_method == "MB") {
      upperbound <- 1 / (2 * pi - 1) * q^2 / N
    }

    # Computing upper-bound of the PFER using approach proposed by SS
    if (PFER_method == "SS") {
      cutoff <- pi
      B <- ceiling(K / 2)
      theta <- q / N
      if (cutoff <= 3 / 4) {
        tmp <- 2 * (2 * cutoff - 1 - 1 / (2 * B))
      } else {
        tmp <- (1 + 1 / B) / (4 * (1 - cutoff + 1 / (2 * B)))
      }
      upperbound <- q^2 / N / tmp

      # Setting to Inf if "out of bounds"
      if ((cutoff < 1 / 2 + min(theta^2, 1 / (2 * B) + 3 / 4 * theta^2)) | (cutoff > 1)) {
        upperbound <- Inf
      }
    }
  } else {
    upperbound <- Inf
  }

  # Re-formatting the upperbound
  if (is.na(upperbound)) {
    upperbound <- Inf
  }

  return(upperbound)
}


#' False Discovery Proportion
#'
#' Computes the (upper-bound) of the FDP as a ratio of the (upper-bound)
#' of the PFER over the number of stably selected features.
#' In stability selection, the FDP corresponds to the expected proportion of stably selected
#' features that are not relevant to the outcome
#' (i.e. proportion of False Positives among stably selected features).
#'
#' @param PFER Per Family Error Rate.
#' @param stab_iter matrix or vector of selection proportions.
#' @param pi vector of thresholds in selection proportions.
#'
#' @return the (upper-bound) of the FDP.
#'
#' @family stability metric functions
#'
#' # Simulating set of selection proportions
#' selprop=round(runif(n=20), digits=2)
#'
#' # Computing the FDP with a threshold of 0.8
#' fdp=FDP(PFER=3, stab_iter=selprop, pi=0.8)
#'
#' @export
FDP <- function(PFER, stab_iter, pi) {
  # Preparing objects
  if (is.matrix(stab_iter)) {
    stab_iter <- stab_iter[upper.tri(stab_iter)]
  }

  # Computing the number of stable edges
  S <- sum(stab_iter >= pi, na.rm = TRUE)

  # Computing the proportion of false discoveries among discoveries (False Discovery Proportion)
  if (S != 0) {
    FDP <- PFER / S
  } else {
    FDP <- 0
  }

  return(FDP)
}

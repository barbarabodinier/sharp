#' Per Family Error Rate
#'
#' Computes the Per Family Error Rate upper-bound of a stability selection model
#' using the methods proposed by Meinshausen and Bühlmann (2010) or Shah and
#' Samworth (2013). In stability selection, the PFER corresponds to the expected
#' number of stably selected features that are not relevant to the outcome (i.e.
#' False Positives).
#'
#' @inheritParams VariableSelection
#' @param q average number of features selected by the underlying algorithm.
#' @param N total number of features.
#' @param pi threshold in selection proportions.
#'
#' @return The estimated upper-bound in PFER.
#'
#' @references \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{stabilityselectionSS}{sharp}
#'
#' @family stability metric functions
#'
#' @examples
#' # Computing PFER for 10/50 selected features and threshold of 0.8
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
#' Computes the False Discovery Proportion (upper-bound) as a ratio of the PFER
#' (upper-bound) over the number of stably selected features. In stability
#' selection, the FDP corresponds to the expected proportion of stably selected
#' features that are not relevant to the outcome (i.e. proportion of False
#' Positives among stably selected features).
#'
#' @param selprop matrix or vector of selection proportions.
#' @param PFER Per Family Error Rate.
#' @param pi threshold in selection proportions.
#'
#' @return The estimated upper-bound in FDP.
#'
#' @family stability metric functions
#'
#' @examples
#' # Simulating set of selection proportions
#' selprop <- round(runif(n = 20), digits = 2)
#'
#' # Computing the FDP with a threshold of 0.8
#' fdp <- FDP(PFER = 3, selprop = selprop, pi = 0.8)
#' @export
FDP <- function(selprop, PFER, pi) {
  # Preparing objects
  if (is.matrix(selprop)) {
    selprop <- selprop[upper.tri(selprop)]
  }

  # Computing the number of stable edges
  S <- sum(selprop >= pi, na.rm = TRUE)

  # Computing the proportion of false discoveries among discoveries (False Discovery Proportion)
  if (S != 0) {
    FDP <- PFER / S
  } else {
    FDP <- 0
  }

  return(FDP)
}

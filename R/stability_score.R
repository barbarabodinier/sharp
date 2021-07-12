#' Stability score
#'
#' Computes the stability score from selection proportions of models with a
#' given parameter controlling the sparsity and for different thresholds in
#' selection proportions. The score measures how unlikely it is that the
#' selection procedure is uniform (i.e. uninformative) for a given combination
#' of parameters.
#'
#' @inheritParams StabilityMetrics
#'
#' @return A vector of stability scores obtained with the different thresholds
#'   in selection proportions.
#'
#' @family stability metric functions
#'
#' @references \insertRef{ourstabilityselection}{focus}
#'
#' @examples
#' # Simulating set of selection proportions
#' selprop <- round(runif(n = 20), digits = 2)
#'
#' # Computing stability scores for different thresholds
#' score <- StabilityScore(selprop, pi_list = c(0.6, 0.7, 0.8), K = 100)
#' @export
StabilityScore <- function(selprop, pi_list = seq(0.6, 0.9, by = 0.01), K, n_cat = 3, group = NULL) {
  # Preparing objects
  if (is.matrix(selprop)) {
    selprop <- selprop[upper.tri(selprop)]
  }

  # Using group penalisation (extracting one per group)
  if (!is.null(group)) {
    selprop <- selprop[cumsum(group)]
  }

  # Computing the number of features (edges/variables)
  N <- length(selprop)

  # Computing the average number of selected features
  q <- round(sum(selprop))

  # Loop over the values of pi
  score <- rep(NA, length(pi_list))
  for (i in 1:length(pi_list)) {
    pi <- pi_list[i]

    # Computing the probabilities of being stable-in, stable-out or unstable under the null (uniform selection)
    p_vect <- BinomialProbabilities(q, N, pi, K, n_cat = n_cat)

    # Computing the log-likelihood
    if (any(is.na(p_vect))) {
      # Returning NA if not possible to compute (e.g. negative number of unstable features, as with pi<=0.5)
      l <- NA
    } else {
      if (n_cat == 2) {
        S_0 <- sum(selprop < pi) # Number of not stably selected features
        S_1 <- sum(selprop >= pi) # Number of stably selected features

        # Checking consistency
        if (S_0 + S_1 != N) {
          stop(paste0("Inconsistency in number of edges \n S_0+S_1=", S_0 + S_1, " instead of ", N))
        }

        # Log-likelihood
        l <- S_0 * p_vect$p_0 + S_1 * p_vect$p_1
      }

      if (n_cat == 3) {
        S_0 <- sum(selprop <= (1 - pi)) # Number of stable-out features
        S_1 <- sum(selprop >= pi) # Number of stable-in features
        U <- sum((selprop < pi) & (selprop > (1 - pi))) # Number of unstable features

        # Checking consistency
        if (S_0 + S_1 + U != N) {
          stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0 + S_1 + U, " instead of ", N))
        }

        # Log-likelihood
        l <- S_0 * p_vect$p_1 + U * p_vect$p_2 + S_1 * p_vect$p_3
      }

      # Re-formatting if infinite
      if (is.infinite(l)) {
        l <- NA
      }
    }

    # Getting the stability score
    score[i] <- -l
  }

  return(score)
}


#' Binomial probabilities for stability score
#'
#' Computes the probabilities of observing each category of selection
#' proportions under the assumption of a uniform selection procedure.
#'
#' @inheritParams StabilityScore
#' @param q average number of features selected by the underlying algorithm.
#' @param N total number of features.
#' @param pi threshold in selection proportions. If n_cat=3, these values must
#'   be >0.5 and <1. If n_cat=2, these values must be >0 and <1.
#'
#' @return A list of probabilities for each of the 2 or 3 categories of
#'   selection proportions.
#'
#' @export
BinomialProbabilities <- function(q, N, pi, K, n_cat = 3) {
  if (n_cat == 2) {
    # Definition of the threshold in selection counts
    thr <- round(K * pi) # Threshold above (>=) which the feature is stably selected

    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_0 <- stats::pbinom(thr - 1, size = K, prob = q / N, log.p = TRUE) # proportion < pi

    # Probability of observing a selection count above thr_up under the null
    p_1 <- stats::pbinom(thr - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE) # proportion >= pi

    # Checking consistency between the three computed probabilities (should sum to 1)
    if (abs(exp(p_0) + exp(p_1) - 1) > 1e-3) {
      message(paste("N:", N))
      message(paste("q:", q))
      message(paste("K:", K))
      message(paste("pi:", pi))
      stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_0+p_1=", exp(p_0) + exp(p_1)))
    }

    # Output the two probabilities under the assumption of uniform selection procedure
    return(list(p_0 = p_0, p_1 = p_1))
  }

  if (n_cat == 3) {
    # Definition of the two thresholds in selection counts
    thr_down <- round(K * (1 - pi)) # Threshold below (<=) which the feature is stable-out
    thr_up <- round(K * pi) # Threshold above (>=) which the feature is stable-in

    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_1 <- stats::pbinom(thr_down, size = K, prob = q / N, log.p = TRUE) # proportion <= (1-pi)

    # Probability of observing a selection count between thr_down and thr_up under the null
    if ((thr_down) >= (thr_up - 1)) {
      # Not possible to compute (i.e. negative number of unstable features)
      p_2 <- NA
    } else {
      # Using cumulative probabilities
      p_2 <- log(stats::pbinom(thr_up - 1, size = K, prob = q / N) - stats::pbinom(thr_down, size = K, prob = q / N)) # 1-pi < proportion < pi

      # Using sum of probabilities (should not be necessary)
      if (is.infinite(p_2) | is.na(p_2)) {
        p_2 <- 0
        for (i in seq(thr_down + 1, thr_up - 1)) {
          p_2 <- p_2 + stats::dbinom(i, size = K, prob = q / N)
        }
        p_2 <- log(p_2)
      }
    }

    # Probability of observing a selection count above thr_up under the null
    p_3 <- stats::pbinom(thr_up - 1, size = K, prob = q / N, lower.tail = FALSE, log.p = TRUE) # proportion >= pi

    # Checking consistency between the three computed probabilities (should sum to 1)
    if (!is.na(p_2)) {
      if (abs(exp(p_1) + exp(p_2) + exp(p_3) - 1) > 1e-3) {
        message(paste("N:", N))
        message(paste("q:", q))
        message(paste("K:", K))
        message(paste("pi:", pi))
        stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1) + exp(p_2) + exp(p_3)))
      }
    }

    # Output the three probabilities under the assumption of uniform selection procedure
    return(list(p_1 = p_1, p_2 = p_2, p_3 = p_3))
  }
}

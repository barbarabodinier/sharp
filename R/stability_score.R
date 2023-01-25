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
#' @details The stability score is derived from the likelihood under the
#'   assumption of uniform (uninformative) selection.
#'
#'   We classify the features into three categories: the stably selected ones
#'   (that have selection proportions \eqn{\ge \pi}), the stably excluded ones
#'   (selection proportion \eqn{\le 1-\pi}), and the unstable ones (selection
#'   proportions between \eqn{1-\pi} and \eqn{\pi}).
#'
#'   Under the hypothesis of equiprobability of selection (instability), the
#'   likelihood of observing stably selected, stably excluded and unstable
#'   features can be expressed as:
#'
#'   \eqn{L_{\lambda, \pi} = \prod_{j=1}^N [ ( 1 - F( K \pi - 1 ) )^{1_{H_{\lambda} (j) \ge K \pi}}
#'   \times ( F( K \pi - 1 ) - F( K ( 1 - \pi ) )^{1_{ (1-\pi) K < H_{\lambda} (j) < K \pi }}
#'   \times F( K ( 1 - \pi ) )^{1_{ H_{\lambda} (j) \le K (1-\pi) }} ]}
#'
#'   where \eqn{H_{\lambda} (j)} is the selection count of feature \eqn{j} and
#'   \eqn{F(x)} is the cumulative probability function of the binomial
#'   distribution with parameters \eqn{K} and the average proportion of selected
#'   features over resampling iterations.
#'
#'   The stability score is computed as the minus log-transformed likelihood
#'   under the assumption of uniform selection:
#'
#'   \eqn{S_{\lambda, \pi} = -log(L_{\lambda, \pi})}
#'
#'   Alternatively, the stability score can be computed by considering only two
#'   sets of features: stably selected (selection proportions \eqn{\ge \pi}) or
#'   not (selection proportions \eqn{< \pi}). This can be done using
#'   \code{n_cat=2}.
#'
#' @return A vector of stability scores obtained with the different thresholds
#'   in selection proportions.
#'
#' @family stability metric functions
#'
#' @references \insertRef{ourstabilityselection}{sharp}
#'
#' @examples
#' # Simulating set of selection proportions
#' set.seed(1)
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
  N <- sum(!is.na(selprop))

  # Computing the average number of selected features
  q <- round(sum(selprop, na.rm = TRUE))

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
        S_0 <- sum(selprop < pi, na.rm = TRUE) # Number of not stably selected features
        S_1 <- sum(selprop >= pi, na.rm = TRUE) # Number of stably selected features

        # Checking consistency
        if (S_0 + S_1 != N) {
          stop(paste0("Inconsistency in number of edges \n S_0+S_1=", S_0 + S_1, " instead of ", N))
        }

        # Log-likelihood
        l <- S_0 * p_vect$p_0 + S_1 * p_vect$p_1
      }

      if (n_cat == 3) {
        S_0 <- sum(selprop <= (1 - pi), na.rm = TRUE) # Number of stable-out features
        S_1 <- sum(selprop >= pi, na.rm = TRUE) # Number of stable-in features
        U <- sum((selprop < pi) & (selprop > (1 - pi)), na.rm = TRUE) # Number of unstable features

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
#' @keywords internal
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


#' Consensus score
#'
#' Computes the consensus score from the consensus matrix. The score measures
#' how unlikely it is that the clustering procedure is uniform (i.e.
#' uninformative) for a given combination of parameters.
#'
#' @param coprop consensus matrix obtained with \code{nc} clusters across
#'   \code{K} subsampling iterations.
#' @param nc number of clusters.
#' @param K number of subsampling iterations.
#' @param linkage character string indicating the type of linkage used in
#'   hierarchical clustering to define the stable clusters. Possible values
#'   include \code{"complete"}, \code{"single"} and \code{"average"} (see
#'   argument \code{"method"} in \code{\link[stats]{hclust}} for a full list).
#'
#' @details Let \eqn{\Gamma(\lambda, G)} be the consensus matrix. We introduce
#'   the matrix \eqn{H(\lambda, G)} of co-membership count corrected for the
#'   subsampling procedure, defined as the integer part of \eqn{K
#'   \Gamma(\lambda, G)}.
#'
#'   Under the hypothesis of equiprobability of co-membership (instability), we
#'   assume that the co-membership counts follow the same binomial distribution
#'   for all pairs of items.
#'
#'   Given stable clusters \eqn{Z} and hyper-parameters \eqn{(\lambda, G)},
#'   clustering stability is measured as the probability \eqn{p_{\lambda,
#'   G}(H|z)} of observing co-membership counts in \eqn{H} that are at least as
#'   high within clusters and at least as low between clusters under
#'   equiprobability:
#'
#'   \eqn{p_{\lambda, G}(H|Z) = \prod_{i < j} F(H_{ij})^{1_{Z_i=/=Z_j}}
#'   \times (1 - F(H_{ij}))^{1_{Z_i=Z_j}}}
#'
#'   where \eqn{F(x)} is the cumulative probability function of the binomial
#'   distribution with parameters \eqn{K} and probability \eqn{\gamma = N_c/N}
#'   with \eqn{N_c} the number of stable co-members and \eqn{N} the number of
#'   item pairs.
#'
#'   This probability is minimised at \eqn{H^s}, which corresponding to the most
#'   stable clustering and is defined as:
#'
#'   \eqn{H^s_{ij} = K \times 1_{Z_i = Z_j}}
#'
#'   The consensus score is calculated as the following standardised
#'   probability:
#'
#'   \eqn{S_c(\lambda, G) = (p_{\lambda, G}(H (\lambda, G)|Z)) / (p_{\lambda,
#'   G}(H^s|Z))}
#'
#'   The consensus score increases with clustering stability.
#'
#' @return A consensus score between 0 and 1.
#'
#' @family stability metric functions
#'
#' @examples
#'
#' # Data simulation
#' set.seed(2)
#' simul <- SimulateClustering(
#'   n = c(30, 30, 30),
#'   nu_xc = 1
#' )
#' plot(simul)
#'
#' # Consensus clustering
#' stab <- Clustering(
#'   xdata = simul$data
#' )
#' stab$Sc[3]
#'
#' # Calculating the consensus score
#' ConsensusScore(
#'   coprop = stab$coprop[, , 3],
#'   nc = stab$nc[3]
#' )
#'
#' @export
ConsensusScore <- function(coprop, nc, K = 100, linkage = "complete") {
  # Clustering on the consensus matrix
  sh_clust <- stats::hclust(stats::as.dist(1 - coprop), method = linkage)

  # Identifying stable clusters
  theta <- stats::cutree(sh_clust, k = nc)

  H <- round(coprop * K)
  H_within <- H * CoMembership(theta)
  N_c <- sum(K * CoMembership(theta)[upper.tri(CoMembership(theta))])

  loglik <- stats::dhyper(
    x = sum(H_within[upper.tri(H_within)]),
    m = N_c,
    n = sum(K * (1 - CoMembership(theta))[upper.tri(CoMembership(theta))]),
    k = sum(H[upper.tri(H)]), log = TRUE
  )

  # max_within=min(sum(H[upper.tri(H)]), N_c)
  # loglik_best=c(loglik_best, dhyper(x=max_within,
  #                                   m = N_c,
  #                                   n = sum(K*(1-CoMembership(theta))[upper.tri(CoMembership(theta))]),
  #                                   k = sum(H[upper.tri(H)]), log = TRUE))

  return(-loglik)
}

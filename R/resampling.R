#' Resampling observations
#'
#' Generates a vector of resampled observation IDs.
#'
#' @param data vector or matrix of data. In regression, this should be the
#'   outcome data.
#' @param family type of regression model. This argument is defined as in the
#'   \code{\link[glmnet]{glmnet}} function from the glmnet package. Possible values
#'   include "gaussian" (linear regression), "binomial" (logistic regression),
#'   "multinomial" (multinomial regression), and "cox" (survival analysis). This
#'   argument is only used with implementation="glmnet", or with functions using
#'   the family argument in the same way.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param resampling resampling approach. Possible values are: "subsampling" for
#'   sampling without replacement of a proportion tau of the observations, or
#'   "bootstrap" for sampling with replacement generating a resampled dataset
#'   with as many observations as in the full sample. Alternatively, this
#'   argument can be a character string indicating the name of a function to use
#'   for resampling. This function must use arguments called "data" and "tau"
#'   and return IDs of observations to be included in the resampled dataset (see
#'   example below).
#' @param ... additional parameters passed to the functions provided in
#'   "resampling".
#'
#' @return a vector of resampled IDs.
#'
#' @details With categorical outcomes (i.e. "family" argument is set to
#'   "binomial", "multinomial" or "cox"), the resampling is done such that the
#'   proportion of observations from each of the categories is representative of
#'   that of the full sample.
#'
#' @examples
#' \dontrun{
#'
#' ## Linear regression framework
#' # Data simulation
#' simul <- SimulateRegression()
#'
#' # Subsampling
#' ids <- Resample(data = simul$Y, family = "gaussian")
#' sum(duplicated(ids))
#'
#' # Bootstrapping
#' ids <- Resample(data = simul$Y, family = "gaussian", resampling = "bootstrap")
#' sum(duplicated(ids))
#'
#' ## Logistic regression framework
#' # Data simulation
#' simul <- SimulateRegression(family = "binomial")
#'
#' # Subsampling
#' ids <- Resample(data = simul$Y, family = "binomial")
#' sum(duplicated(ids))
#' prop.table(table(simul$Y))
#' prop.table(table(simul$Y[ids]))
#'
#' # Data simulation for a binary confounder
#' conf <- ifelse(runif(n = 100) > 0.5, yes = 1, no = 0)
#'
#' # User-defined resampling function
#' BalancedResampling <- function(data, tau, Z, ...) {
#'   s <- NULL
#'   for (z in unique(Z)) {
#'     s <- c(s, sample(which((data == "0") & (Z == z)), size = tau * sum((data == "0") & (Z == z))))
#'     s <- c(s, sample(which((data == "1") & (Z == z)), size = tau * sum((data == "1") & (Z == z))))
#'   }
#'   return(s)
#' }
#'
#' # Resampling keeping proportions by Y and Z
#' ids <- Resample(data = simul$Y, family = "binomial", resampling = "BalancedResampling", Z = conf)
#' prop.table(table(simul$Y, conf))
#' prop.table(table(simul$Y[ids], conf[ids]))
#'
#' # User-defined resampling for stability selection
#' stab <- VariableSelection(
#'   xdata = simul$X, ydata = simul$Y, family = "binomial",
#'   resampling = "BalancedResampling", Z = conf
#' )
#' }
#' @export
Resample <- function(data, family = NULL, tau = 0.5, resampling = "subsampling", ...) {
  # Preparing the data
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }

  if (!resampling %in% c("subsampling", "bootstrap")) {
    s <- do.call(get(resampling), args = list(data = data, tau = tau, ...))
  } else {
    # Using or not replacement in resampling
    replacement <- ifelse(resampling == "subsampling", yes = FALSE, no = TRUE)

    # Definition of the size of sub/bootstrap sample
    if (replacement) {
      tau <- 1
    }

    # Resampling procedure
    if (!is.null(family)) {
      # Resampling for regression models
      if (family %in% c("gaussian", "poisson", "mgaussian")) {
        s <- sample(nrow(data), size = tau * nrow(data), replace = replacement)
      }
      if (family == "binomial") {
        data <- cbind(apply(data, 1, sum)) # to ensure balanced classes for PLS-DA
        s <- NULL
        for (mycat in levels(factor(data))) {
          scat <- sample(which(data == mycat), size = tau * sum(data == mycat), replace = replacement)
          s <- c(s, scat)
        }
      }
      if (family == "multinomial") {
        s <- NULL
        for (mycat in levels(factor(data))) {
          scat <- sample(which(data == mycat), size = tau * sum(data == mycat), replace = replacement)
          s <- c(s, scat)
        }
      }
      if (family == "cox") {
        s0 <- sample(which(data[, 2] == "0"), size = tau * sum(data[, 2] == "0"), replace = replacement)
        s1 <- sample(which(data[, 2] == "1"), size = tau * sum(data[, 2] == "1"), replace = replacement)
        s <- c(s0, s1)
      }
    } else {
      # Resampling for network models
      s <- sample(1:nrow(data), size = tau * nrow(data), replace = replacement)
    }
  }
  return(s)
}

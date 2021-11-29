#' Resampling observations
#'
#' Generates a vector of resampled observation IDs.
#'
#' @inheritParams VariableSelection
#' @param data vector or matrix of data. In regression, this should be the
#'   outcome data.
#' @param ... additional parameters passed to the function provided in
#'   \code{resampling}.
#'
#' @return A vector of resampled IDs.
#'
#' @details With categorical outcomes (i.e. "family" argument is set to
#'   "binomial", "multinomial" or "cox"), the resampling is done such that the
#'   proportion of observations from each of the categories is representative of
#'   that of the full sample.
#'
#' @examples
#'
#' ## Linear regression framework
#' # Data simulation
#' simul <- SimulateRegression()
#'
#' # Subsampling
#' ids <- Resample(data = simul$ydata, family = "gaussian")
#' sum(duplicated(ids))
#'
#' # Bootstrapping
#' ids <- Resample(data = simul$ydata, family = "gaussian", resampling = "bootstrap")
#' sum(duplicated(ids))
#'
#' ## Logistic regression framework
#' # Data simulation
#' simul <- SimulateRegression(family = "binomial")
#'
#' # Subsampling
#' ids <- Resample(data = simul$ydata, family = "binomial")
#' sum(duplicated(ids))
#' prop.table(table(simul$ydata))
#' prop.table(table(simul$ydata[ids]))
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
#' ids <- Resample(data = simul$ydata, family = "binomial", resampling = BalancedResampling, Z = conf)
#' prop.table(table(simul$ydata, conf))
#' prop.table(table(simul$ydata[ids], conf[ids]))
#'
#' # User-defined resampling for stability selection
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata, family = "binomial",
#'   resampling = BalancedResampling, Z = conf
#' )
#' @export
Resample <- function(data, family = NULL, tau = 0.5, resampling = "subsampling", ...) {
  # Preparing the data
  if (is.factor(data)) {
    data <- as.character(factor(data, levels = levels(data), labels = seq(1, length(levels(data))) - 1))
  }
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.null(family)){
    if (family == "multinomial") {
      if (is.matrix(data)) {
        data <- DummyToCategories(x = data, verbose = FALSE)
      }
    }
  }
  
  # if (!resampling %in% c("subsampling", "bootstrap")) {
  if (is.function(resampling)) {
    # s <- do.call(get(resampling), args = list(data = data, tau = tau, ...))
    s <- do.call(resampling, args = list(data = data, tau = tau, ...))
  } else {
    if (!resampling %in% c("subsampling", "bootstrap")) {
      stop("Invalid input for argument 'resampling'. It must be a function or a character string: 'subsampling' or 'bootstrap'.")
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
          if (ncol(data) > 1) {
            data <- cbind(apply(data, 1, sum)) # to ensure balanced classes for PLS-DA
          }
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
  }
  return(s)
}


#' Splitting observations into folds
#'
#' Generates a list of \code{n_folds} non-overlapping sets of observation IDs
#' (folds).
#'
#' @inheritParams Resample
#' @param n_folds number of folds.
#'
#' @return A list of length \code{n_folds} with sets of non-overlapping
#'   observation IDs.
#'
#' @details With categorical outcomes (i.e. "family" argument is set to
#'   "binomial", "multinomial" or "cox"), the split is done such that the
#'   proportion of observations from each of the categories in each of the folds
#'   is representative of that of the full sample.
#'
#' @examples
#'
#' # Splitting into 5 folds
#' simul <- SimulateRegression()
#' ids <- Folds(data = simul$ydata)
#' lapply(ids, length)
#'
#' # Balanced folds with respect to a binary variable
#' simul <- SimulateRegression(family = "binomial")
#' ids <- Folds(data = simul$ydata, family = "binomial")
#' lapply(ids, FUN = function(x) {
#'   table(simul$ydata[x, ])
#' })
#' @export
Folds <- function(data, family = NULL, n_folds = 5) {
  # Re-formatting the inputs
  if (is.vector(data)) {
    data <- cbind(data)
  }
  rownames(data) <- paste0("obs", 1:nrow(data))
  
  # Storing total number of observations
  n <- nrow(data)
  
  # Creating balanced folds
  folds_ids <- list()
  for (k in 1:n_folds) {
    ids <- rownames(data)[Resample(data = data, family = family, tau = 1 / n_folds * n / nrow(data))]
    folds_ids <- c(folds_ids, list(ids))
    data <- data[-which(rownames(data) %in% ids), , drop = FALSE]
  }
  
  # Returning row ids
  folds_ids <- lapply(folds_ids, FUN = function(x) {
    as.numeric(gsub("obs", "", x))
  })
  
  return(folds_ids)
}

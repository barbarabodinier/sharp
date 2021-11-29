#' True and False Positive Rates
#'
#' Computes the True and False Positive Rates by comparing the true (observed) and
#' predicted status. The predicted status is obtained by applying a threshold on
#' the predicted scores.
#'
#' @param observed observed binary status.
#' @param predicted predicted score.
#' @param thr threshold for predicted probabilities.
#'
#' @return True and False Positive Rates (TPR and FPR, respectively).
#'
#' @keywords internal
Rates <- function(observed, predicted, thr) {
  contingency <- table(
    factor(predicted > thr, levels = c(FALSE, TRUE)),
    factor(observed, levels = c(0, 1))
  )

  TP <- contingency[2, 2]
  P <- sum(contingency[, 2])
  TPR <- TP / P

  FP <- contingency[2, 1]
  N <- sum(contingency[, 1])
  FPR <- FP / N

  return(list(TPR = TPR, FPR = FPR))
}


#' Receiver Operating Characteristic (ROC)
#'
#' Computes the True and False Positive Rates (TPR and FPR, respectively) and
#' Area Under the Curve (AUC) by comparing the true (observed) and predicted
#' status using a range of thresholds on the predicted score.
#'
#' @param predicted numeric predicted scores.
#' @param observed factor encoding the observed binary status.
#' @param n_thr number of thresholds to use to construct the ROC curve. For
#'   faster computations on large data, values below \code{length(x)-1} can be
#'   used.
#'
#' @return A list with: \item{TPR}{True Positive Rate.} \item{FPR}{False
#'   Positive Rate.} \item{AUC}{Area Under the Curve.}
#'
#' @family prediction performance functions
#'
#' @examples
#' \dontrun{
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 20, family = "binomial")
#'
#' # Balanced training/test split
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' x2 <- simul$xdata[-ids_train, , drop = FALSE]
#' y2 <- simul$ydata[-ids_train, , drop = FALSE]
#' ids_recalib <- Resample(
#'   data = y2,
#'   tau = 0.5, family = "binomial"
#' )
#' xrecalib <- x2[ids_recalib, , drop = FALSE]
#' yrecalib <- y2[ids_recalib, , drop = FALSE]
#' xtest <- x2[-ids_recalib, ]
#' ytest <- y2[-ids_recalib, ]
#'
#' # Stability selection and recalibration
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#' recalibrated <- Recalibrate(xdata = xrecalib, ydata = yrecalib, stability = stab)
#'
#' # ROC analysis
#' predicted <- predict(recalibrated, newdata = as.data.frame(xtest))
#' roc <- ROC(predicted = predicted, observed = ytest)
#' PlotROC(roc)
#' }
#'
#' @export
ROC <- function(predicted, observed, n_thr = NULL) {
  # Checking the inputs
  predicted <- as.numeric(predicted)
  if (is.factor(observed)) {
    observed <- factor(observed, levels = levels(observed), labels = c(0, 1))
  } else {
    observed <- factor(observed, levels = sort(unique(observed)), labels = c(0, 1))
  }

  # Defining the thresholds
  breaks <- sort(unique(predicted), decreasing = FALSE)
  breaks <- breaks[-1]
  if (!is.null(n_thr)) {
    if (length(breaks) > n_thr) {
      breaks <- breaks[floor(seq(1, length(breaks), length.out = n_thr))]
    } else {
      breaks <- seq(min(breaks), max(breaks), length.out = n_thr)
    }
  }

  # Computing
  TPR <- FPR <- rep(NA, length(breaks) + 2)
  for (k in 1:length(breaks)) {
    out <- Rates(observed = observed, predicted = predicted, thr = breaks[k])
    TPR[k + 1] <- out$TPR
    FPR[k + 1] <- out$FPR
  }
  TPR[1] <- FPR[1] <- 1
  TPR[length(TPR)] <- FPR[length(FPR)] <- 0

  # Computing the AUC
  tmp <- apply(rbind(TPR[-1], TPR[-length(TPR)]), 2, mean)
  AUC <- abs(sum(diff(FPR) * tmp))

  return(list(FPR = rbind(FPR), TPR = rbind(TPR), AUC = AUC))
}


#' Regression model recalibration
#'
#' Recalibrates the regression using an un-penalised model with stably selected
#' variables as predictors. Variables in \code{xdata} not evaluated in the
#' stability selection model will automatically be included as predictors.
#'
#' @inheritParams VariableSelection
#' @param stability output of \code{\link{VariableSelection}}. With
#'   \code{stability=NULL} (the default), a model including all variables in
#'   \code{xdata} as predictors is fitted. Argument \code{family} must be
#'   provided in this case.
#' @param family type of regression model. Possible values include
#'   \code{"gaussian"} (linear regression), \code{"binomial"} (logistic
#'   regression), \code{"multinomial"} (multinomial regression), and
#'   \code{"cox"} (survival analysis). This argument must be consistent with
#'   input \code{stability}, if provided.
#'
#' @return The output as obtained from: \item{\code{\link[stats]{lm}}}{for
#'   linear regression (\code{"gaussian"} family).}
#'   \item{\code{\link[stats]{glm}}}{for logistic regression (\code{"binomial"}
#'   family).} \item{\code{\link[nnet]{multinom}}}{for multinomial regression
#'   (\code{"multinomial"} family).}
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \dontrun{
#'
#' ## Linear regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#'
#' # Data split
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "gaussian"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' xrecalib <- simul$xdata[-ids_train, , drop = FALSE]
#' yrecalib <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "gaussian")
#' print(SelectedVariables(stab))
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = xrecalib, ydata = yrecalib,
#'   stability = stab
#' )
#' recalibrated$coefficients # recalibrated coefficients
#' head(recalibrated$fitted.values) # recalibrated predicted values
#'
#' # Fitting the full model (including all possible predictors)
#' recalibrated <- Recalibrate(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian"
#' )
#' recalibrated$coefficients # recalibrated coefficients
#'
#'
#' ## Cox regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "binomial")
#' ydata <- cbind(
#'   time = runif(nrow(simul$ydata), min = 100, max = 2000),
#'   case = simul$ydata[, 1]
#' ) # including dummy time to event
#'
#' # Data split
#' ids_train <- Resample(
#'   data = ydata,
#'   tau = 0.5, family = "cox"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- ydata[ids_train, , drop = FALSE]
#' xrecalib <- simul$xdata[-ids_train, , drop = FALSE]
#' yrecalib <- ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "cox")
#' print(SelectedVariables(stab))
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = xrecalib, ydata = yrecalib,
#'   stability = stab
#' )
#' recalibrated$coefficients # recalibrated coefficients
#' head(recalibrated$linear.predictors) # recalibrated scores
#'
#'
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#'
#' # Data split
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' xrecalib <- simul$xdata[-ids_train, , drop = FALSE]
#' yrecalib <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = xrecalib, ydata = yrecalib,
#'   stability = stab
#' )
#' recalibrated$coefficients # recalibrated coefficients
#' head(recalibrated$fitted.values) # recalibrated predicted probabilities
#'
#'
#' ## Multinomial regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 15, family = "multinomial")
#'
#' # Data split
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "multinomial"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' xrecalib <- simul$xdata[-ids_train, , drop = FALSE]
#' yrecalib <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   family = "multinomial"
#' )
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = xrecalib, ydata = yrecalib,
#'   stability = stab
#' )
#' summary(recalibrated) # recalibrated coefficients
#' head(recalibrated$fitted.values) # recalibrated predicted probabilities
#' }
#'
#' @export
Recalibrate <- function(xdata, ydata, stability = NULL, family = NULL) {
  if (!is.null(stability)) {
    if (class(stability) != "variable_selection") {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!stability$methods$family %in% c("gaussian", "cox", "binomial", "multinomial")) {
      stop("This function can only be applied with the following families: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
    }
  }

  # Re-formatting the input
  if (is.vector(xdata)) {
    xdata <- cbind(xdata)
    colnames(xdata) <- "var"
  }

  # Defining predictors for the model
  if (is.null(stability)) {
    selected <- rep(1, ncol(xdata))
    names(selected) <- colnames(xdata)
  } else {
    selected <- SelectedVariables(stability)
  }
  ids <- c(
    names(selected)[which(selected == 1)],
    colnames(xdata)[!colnames(xdata) %in% names(selected)]
  )
  myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", ids, "`"), collapse = " + ")))

  # Recalibration for linear regression
  if (family == "gaussian") {
    mymodel <- stats::lm(myformula, data = as.data.frame(xdata))
  }

  # Recalibration for Cox regression
  if (family == "cox") {
    ydata <- survival::Surv(ydata[, "time"], ydata[, "case"])
    mymodel <- survival::coxph(myformula, data = as.data.frame(xdata))
  }

  # Recalibration for logistic regression
  if (family == "binomial") {
    mymodel <- stats::glm(myformula,
      data = as.data.frame(xdata),
      family = stats::binomial(link = "logit")
    )
  }

  # Recalibration for multinomial regression
  if (family == "multinomial") {
    ydata <- DummyToCategories(ydata)
    mymodel <- nnet::multinom(myformula, data = as.data.frame(xdata))
  }

  return(mymodel)
}


#' Prediction performance in regression
#'
#' Computes the prediction performance for logistic or Cox regression models.
#' This is done by (i) recalibrating the model on a training set including a
#' proportion \code{tau} of the observations, and (ii) evaluating the
#' performance on the remaining observations (test set). For more reliable
#' results, the procedure can be repeated \code{K} times (default \code{K=1}).
#'
#' @inheritParams Recalibrate
#' @param K number of subsampling iterations.
#' @param tau proportion of observations used in the training set.
#' @param seed value of the seed to ensure reproducibility of the results.
#' @param n_thr number of thresholds to use to construct the ROC curve. For
#'   faster computations on large data, values below \code{length(x)-1} can be
#'   used. Only used for logistic regression.
#' @param time numeric indicating the time for which the survival probabilities
#'   are computed. Only used for Cox regression.
#' @param ij_method logical indicating if the analysis should be done for only
#'   one recalibration/test split with variance of the concordance index should
#'   be computed using the infinitesimal jackknife method as implemented in
#'   \code{\link[survival]{concordance}}. With \code{ij_method=TRUE}, the
#'   concordance index and estimated confidence interval at level 0.05 are
#'   reported. With \code{ij_method=FALSE} (the default), the concordance
#'   indices computed for different recalibration/test splits are reported.
#'
#' @details For a fair evaluation of the prediction performance, the data is
#'   split into a training set (including a proportion \code{tau} of the
#'   observations) and test set (remaining observations). The regression model
#'   is fitted on the training set and applied on the test set. Performance
#'   metrics are computed in the test set by comparing predicted and observed
#'   outcomes.
#'
#'   For logistic regression, a Receiver Operating Characteristic (ROC) analysis
#'   is performed: the True and False Positive Rates (TPR and FPR) are computed
#'   for different thresholds in predicted probabilities.
#'
#'   For Cox regression, the Concordance Index (as implemented in
#'   \code{\link[survival]{concordance}}) looking at survival probabilities up
#'   to a specific \code{time} is computed.
#'
#' @return A list with: \item{TPR}{True Positive Rate (for logistic regression
#'   only).} \item{FPR}{False Positive Rate (for logistic regression only).}
#'   \item{AUC}{Area Under the Curve (for logistic regression only).}
#'   \item{concordance}{Concordance index (for Cox regression only).}
#'   \item{lower}{lower bound of the confidence interval at level 0.05 for the
#'   concordance index calculated using the infinitesimal jackknife (for Cox
#'   regression and with \code{ij_method=TRUE}).} \item{upper}{upper bound of
#'   the confidence interval at level 0.05 for the concordance index calculated
#'   using the infinitesimal jackknife (for Cox regression and with
#'   \code{ij_method=TRUE}).}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Recalibrate}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \dontrun{
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 10, family = "binomial")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Evaluation of the performances on recalibrated models (K=1)
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, n_thr = NULL
#' )
#' PlotROC(roc)
#'
#' # Using more recalibration/test splits
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 100
#' )
#' boxplot(roc$AUC, ylab = "AUC")
#' PlotROC(roc)
#'
#' # Comparison with saturated model
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   family = "binomial", K = 100
#' )
#' PlotROC(roc, col = "blue", col_band = "blue", add = TRUE)
#'
#'
#' ## Cox regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 50, family = "binomial")
#' ydata <- cbind(
#'   time = runif(nrow(simul$ydata), min = 100, max = 2000),
#'   case = simul$ydata[, 1]
#' ) # including dummy time to event
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "cox")
#'
#' # Evaluation of the performances on recalibrated models (K=1)
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, ij_method = TRUE
#' )
#' print(perf)
#'
#' # Using more recalibration/test splits
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 10, time = 1000
#' )
#' boxplot(perf$concordance)
#' }
#'
#' @export
ExplanatoryPerformance <- function(xdata, ydata,
                                   stability = NULL, family = NULL,
                                   K = 1, tau = 0.8, seed = 1,
                                   n_thr = NULL,
                                   ij_method = FALSE, time = 1000) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (class(stability) != "variable_selection") {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!stability$methods$family %in% c("cox", "binomial")) {
      stop("This function can only be applied with the following families: 'binomial' or 'cox'.")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
    }
  }

  # Re-formatting input data
  if (is.vector(ydata)) {
    ydata <- cbind(ydata)
  }
  if (ij_method) {
    K <- 1
  }
  if (family == "binomial") {
    metric <- "roc"
  } else {
    metric <- "concordance"
  }
  # if (is.null(n_folds) | (K==1)) {
  #   n_folds <- 1
  # } else {
  #   K <- ceiling(K / n_folds)
  # }
  n_folds <- 1

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Running the subsampling iterations
  iter <- 0
  for (k in 1:K) {
    for (fold_id in 1:n_folds) {
      iter <- iter + 1
      if (n_folds == 1) {
        # Balanced training/test split
        ids_test <- Resample(data = ydata, tau = 1 - tau, family = family)
      } else {
        if (fold_id == 1) {
          ids_folds <- Folds(data = ydata, n_folds = n_folds)
        }
        ids_test <- ids_folds[[fold_id]]
      }
      xtrain <- xdata[-ids_test, , drop = FALSE]
      ytrain <- ydata[-ids_test, , drop = FALSE]
      xtest <- xdata[ids_test, , drop = FALSE]
      ytest <- ydata[ids_test, , drop = FALSE]

      # Recalibration from stability selection model
      recalibrated <- Recalibrate(xdata = xtrain, ydata = ytrain, stability = stability, family = family)

      # Performing ROC analyses
      if (tolower(metric) == "roc") {
        # ROC analysis
        predicted <- stats::predict.glm(recalibrated, newdata = as.data.frame(xtest), type = "response")
        roc <- ROC(predicted = predicted, observed = ytest, n_thr = n_thr)

        # Initialisation of the object
        if (iter == 1) {
          n_thr <- length(roc$FPR) - 2
          FPR <- TPR <- matrix(NA, nrow = K * n_folds, ncol = length(roc$TPR))
          AUC <- rep(NA, K * n_folds)
        }

        # Storing the metrics
        FPR[iter, ] <- roc$FPR
        TPR[iter, ] <- roc$TPR
        AUC[iter] <- roc$AUC
      }

      # Performing concordance analyses
      if (tolower(metric) == "concordance") {
        # Computing the concordance index for given times
        predicted <- stats::predict(recalibrated, newdata = as.data.frame(xtest), type = "lp")
        survobject <- survival::Surv(ytest[, "time"], ytest[, "case"])
        S0 <- summary(survival::survfit(recalibrated), times = time, extend = TRUE)$surv
        S <- S0^exp(predicted)
        cstat <- survival::concordance(survobject ~ S)

        if (ij_method) {
          # Storing the concordance index and confidence interval
          cindex <- cstat$concordance
          lower <- cindex - 1.96 * sqrt(cstat$var)
          upper <- cindex + 1.96 * sqrt(cstat$var)
        } else {
          # Storing concordance index
          if (iter == 1) {
            cindex <- rep(NA, K * n_folds)
          }
          cindex[iter] <- cstat$concordance
        }
      }
    }
  }

  if (tolower(metric) == "roc") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }

  if (tolower(metric) == "concordance") {
    if (ij_method) {
      out <- list(concordance = cindex, lower = lower, upper = upper)
    } else {
      out <- list(concordance = cindex)
    }
  }

  return(out)
}


#' Incremental prediction performance in regression
#'
#' Computes the prediction performance for logistic or Cox regression models
#' where predictors are sequentially added by order of decreasing selection
#' proportion. This function can be used to evaluate the marginal contribution
#' of each of the selected predictors over and above more stable predictors.
#' Performances are evaluated as in \code{\link{ExplanatoryPerformance}}.
#'
#' @inheritParams ExplanatoryPerformance
#' @param n_predictors number of predictors to consider.
#'
#' @return For logistic regression, a list with: \item{FPR}{A list with, for
#'   each of the models (sequentially added predictors), the False Positive
#'   Rates for different thresholds (columns) and different data splits (rows).}
#'   \item{TPR}{A list with, for each of the models (sequentially added
#'   predictors), the True Positive Rates for different thresholds (columns) and
#'   different data splits (rows).} \item{AUC}{A list with, for each of the
#'   models (sequentially added predictors), a vector of Area Under the Curve
#'   (AUC) values obtained with different data splits.} \item{names}{Names of
#'   the predictors by order of inclusion.}
#'
#'   For Cox models, a list with: \item{concordance}{If \code{ij_method=FALSE},
#'   a list with, for each of the models (sequentially added predictors), a
#'   vector of concordance indices obtained with different data splits. If
#'   \code{ij_method=TRUE}, a vector of concordance indices for each of the
#'   models (sequentially added predictors).} \item{lower}{A vector of the lower
#'   bound of the confidence interval at level 0.05 for concordance indices for
#'   each of the models (sequentially added predictors). Only returned if
#'   \code{ij_method=TRUE}.} \item{upper}{A vector of the upper bound of the
#'   confidence interval at level 0.05 for concordance indices for each of the
#'   models (sequentially added predictors). Only returned if
#'   \code{ij_method=TRUE}.} \item{names}{Names of the predictors by order of
#'   inclusion.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Recalibrate}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \dontrun{
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "binomial")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Evaluating marginal contribution of the predictors
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#' PlotIncremental(perf)
#'
#'
#' ## Cox regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "binomial")
#' ydata <- cbind(
#'   time = runif(nrow(simul$ydata), min = 100, max = 2000),
#'   case = simul$ydata[, 1]
#' ) # including dummy time to event
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "cox")
#'
#' # Marginal contribution
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#' PlotIncremental(perf)
#'
#' # Faster computations on a single data split
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, ij_method = TRUE)
#' PlotIncremental(perf)
#' }
#'
#' @export
Incremental <- function(xdata, ydata,
                        stability = NULL, family = NULL,
                        n_predictors = NULL,
                        K = 100, tau = 0.8, seed = 1,
                        n_thr = NULL,
                        ij_method = FALSE, time = 1000) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (class(stability) != "variable_selection") {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!stability$methods$family %in% c("cox", "binomial")) {
      stop("This function can only be applied with the following families: 'binomial' or 'cox'.")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
    }
  }

  # Defining the number of predictors
  if (is.null(n_predictors)) {
    if (!is.null(stability)) {
      # Stopping at the calibrated model
      n_predictors <- sum(SelectedVariables(stability))

      # Adding the variables that are forced in the model with penalty.factor
      n_predictors <- n_predictors + sum(!colnames(xdata) %in% names(SelectedVariables(stability)))
    } else {
      n_predictors <- ncol(xdata)
    }
    if (n_predictors == 0) {
      n_predictors <- 10
    }
  }
  n_predictors <- min(n_predictors, ncol(xdata))

  # Defining the order of inclusion in the model
  if (is.null(stability)) {
    myorder <- colnames(xdata) # order of the columns is used
  } else {
    # Including variables by order of decreasing selection proportions
    myorder <- names(SelectionProportions(stability))[sort.list(SelectionProportions(stability), decreasing = TRUE)]

    # Including the variables that are forced in the model first by order of columns in the data
    myorder <- c(colnames(xdata)[!colnames(xdata) %in% names(SelectedVariables(stability))], myorder)
  }

  # Initialisation of the objects
  if (family == "binomial") {
    TPR <- FPR <- AUC <- list()
  }
  if (family == "cox") {
    if (ij_method) {
      concordance <- lower <- upper <- NULL
    } else {
      concordance <- list()
    }
  }

  for (k in 1:n_predictors) {
    perf <- ExplanatoryPerformance(
      xdata = xdata[, myorder[1:k], drop = FALSE],
      ydata = ydata,
      stability = NULL,
      family = family,
      K = K, tau = tau, seed = seed,
      n_thr = n_thr,
      ij_method = ij_method, time = time
    )
    if (family == "binomial") {
      FPR <- c(FPR, list(perf$FPR))
      TPR <- c(TPR, list(perf$TPR))
      AUC <- c(AUC, list(perf$AUC))
    }
    if (family == "cox") {
      if (ij_method) {
        concordance <- c(concordance, perf$concordance)
        lower <- c(lower, perf$lower)
        upper <- c(upper, perf$upper)
      } else {
        concordance <- c(concordance, list(perf$concordance))
      }
    }
  }

  # Preparing the output
  if (family == "binomial") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }
  if (family == "cox") {
    if (ij_method) {
      out <- list(concordance = concordance, lower = lower, upper = upper)
    } else {
      out <- list(concordance = concordance)
    }
  }

  # Adding variable names
  out <- c(out, names = list(myorder[1:n_predictors]))

  return(out)
}


#' Receiver Operating Characteristic (ROC) curve
#'
#' Plots the True Positive Rate (TPR) as a function of the False Positive Rate
#' (FPR) for different thresholds in predicted probabilities. If the results
#' from multiple ROC analyses are provided (e.g. output from
#' \code{\link{ExplanatoryPerformance}} with large \code{K}), the point-wise
#' median is represented and flanked by a transparent band defined by
#' point-wise \code{quantiles}.
#'
#' @param roc output from \code{\link{ROC}} or
#'   \code{\link{ExplanatoryPerformance}}.
#' @param col colour of the point-wise median curve.
#' @param col_band colour of the band defined by point-wise \code{quantiles}.
#' @param alpha level of opacity for the band.
#' @param lwd line width for the point-wise median curve, as in
#'   \code{\link{par}}.
#' @param lty line type for the point-wise median curve, as in
#'   \code{\link{par}}.
#' @param quantiles point-wise quantiles of the performances defining the band.
#' @param add logical indicating if the curve should be added to the current
#'   plot.
#'
#' @return A plot.
#'
#' @family prediction performance functions
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Recalibrate}}
#'
#' @examples
#' \dontrun{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 10, family = "binomial")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Evaluation of the performances on recalibrated models (K=1)
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, n_thr = NULL
#' )
#' PlotROC(roc)
#'
#' # Using more recalibration/test splits
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 100
#' )
#' PlotROC(roc)
#'
#' # Comparison with saturated model
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   family = "binomial", K = 100
#' )
#' PlotROC(roc, col = "blue", col_band = "blue", add = TRUE)
#' }
#'
#' @export
PlotROC <- function(roc,
                    col = "red", col_band = NULL,
                    alpha = 0.5,
                    lwd = 1, lty = 1,
                    quantiles = c(0.05, 0.95),
                    add = FALSE) {
  # Extracting the number of iterations
  niter <- length(roc$AUC)

  # Defining the band colour
  if (is.null(col_band)) {
    col_band <- col
  }

  # Initialising the plot
  if (!add) {
    plot(NULL,
      xlim = c(0, 1), ylim = c(0, 1),
      type = "l", lwd = 2,
      xlab = "False Positive Rate", ylab = "True Positive Rate", las = 1, cex.lab = 1.3,
      panel.first = graphics::abline(0, 1, lty = 3)
    )
  }

  # Defining quantile bands
  if (nrow(roc$FPR) > 1) {
    xseq <- apply(roc$FPR, 2, FUN = function(x) {
      sort(x)[rev(quantiles) * niter]
    })
    yseq <- apply(roc$TPR, 2, FUN = function(x) {
      sort(x)[quantiles * niter]
    })
    graphics::polygon(c(xseq[1, ], rev(xseq[2, ])),
      c(yseq[1, ], rev(yseq[2, ])),
      col = grDevices::adjustcolor(col_band, alpha.f = alpha),
      border = NA
    )
  }

  # Adding the point-wise average
  graphics::lines(apply(roc$FPR, 2, stats::median), apply(roc$TPR, 2, stats::median),
    type = "l", lwd = lwd, lty = lty, col = col
  )
}


#' Visualisation of incremental performance
#'
#' Represents prediction performances upon sequential inclusion of the
#' predictors in a logistic or Cox regression model as produced by
#' \code{\link{Incremental}}.
#'
#' @param perf output from \code{\link{Incremental}}.
#' @param quantiles quantiles defining the lower and upper bounds.
#' @param ylab label of the y-axis.
#' @param pch type of point, as in \code{\link{points}}.
#' @param cex size of point.
#' @param cex.lab size of label on the y-axis.
#' @param sfrac size of the end bars, as in \code{\link[plotrix]{plotCI}}.
#' @param xlas orientation of labels on the x-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylas orientation of labels on the y-axis, as \code{las} in
#'   \code{\link[graphics]{par}}.
#' @param ylim displayed range along the y-axis.
#' @param bty character string indicating if the box around the plot should be
#'   drawn. Possible values include: \code{"o"} (default, the box is drawn), or
#'   \code{"n"} (no box).
#' @param xgrid logical indicating if a vertical grid should be drawn.
#' @param ygrid logical indicating if a horizontal grid should be drawn.
#' @param col vector of point colours.
#' @param col.axis optional vector of label colours. If \code{col.axis=NULL},
#'   the colours provided in argument \code{col} are used.
#' @param xcex.axis size of labels along the x-axis.
#' @param ycex.axis size of labels along the y-axis.
#'
#' @return A plot.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Recalibrate}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \dontrun{
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "binomial")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Evaluating marginal contribution of the predictors
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#'
#' # Basic visualisation
#' PlotIncremental(perf)
#'
#' # Adding grids
#' PlotIncremental(perf, xgrid = TRUE, ygrid = TRUE)
#'
#' # Changing colours
#' PlotIncremental(perf,
#'   bty = "n",
#'   col = colorRampPalette(c("blue", "red"))(length(perf$names))
#' )
#'
#'
#' ## Cox regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "binomial")
#' ydata <- cbind(
#'   time = runif(nrow(simul$ydata), min = 100, max = 2000),
#'   case = simul$ydata[, 1]
#' ) # including dummy time to event
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "binomial"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "cox")
#'
#' # Marginal contribution
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#' PlotIncremental(perf)
#'
#' # Faster computations on a single data split
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, ij_method = TRUE)
#' PlotIncremental(perf)
#' }
#'
#' @export
PlotIncremental <- function(perf, quantiles = c(0.05, 0.95),
                            ylab = "Performance",
                            pch = 18,
                            col = "black", col.axis = NULL,
                            cex = 1, cex.lab = 1.5,
                            xcex.axis = 1, ycex.axis = 1,
                            xlas = 2, ylas = 1,
                            sfrac = 0.005,
                            ylim = NULL, bty = "o",
                            xgrid = FALSE, ygrid = FALSE) {
  # Checking plotrix package is installed
  CheckPackageInstalled("plotrix")

  # Checking the inputs
  quantiles <- sort(quantiles)

  # Re-formatting the inputs
  if (is.null(col.axis)) {
    col.axis <- col
  }

  if ("concordance" %in% names(perf)) {
    if ("lower" %in% names(perf)) {
      x <- perf$concordance
      xlower <- perf$lower
      xupper <- perf$upper
    } else {
      x <- sapply(perf$concordance, mean)
      xlower <- sapply(perf$concordance, stats::quantile, probs = quantiles[1])
      xupper <- sapply(perf$concordance, stats::quantile, probs = quantiles[2])
    }
  } else {
    x <- sapply(perf$AUC, mean)
    xlower <- sapply(perf$AUC, stats::quantile, probs = quantiles[1])
    xupper <- sapply(perf$AUC, stats::quantile, probs = quantiles[2])
  }
  xseq <- 1:length(x)

  # Re-formatting label colours
  if (length(col.axis) < length(x)) {
    col.axis <- rep(col.axis, length(x))
  }

  # Defining the plot range
  if (is.null(ylim)) {
    ylim <- range(c(xlower, xupper, x))
  }

  # Defining horizontal grid
  hseq <- NULL
  if (ygrid) {
    hseq <- grDevices::axisTicks(ylim, log = FALSE)
  }

  # Defining vertical grid
  vseq <- NULL
  if (xgrid) {
    vseq <- xseq
  }

  # Creating the plot
  plot(NULL,
    cex.lab = cex.lab,
    bty = bty,
    xlim = range(xseq), ylim = ylim,
    panel.first = c(
      graphics::abline(h = hseq, col = "grey", lty = 3),
      graphics::abline(v = vseq, col = "grey", lty = 3)
    ),
    xlab = "", xaxt = "n", yaxt = "n", ylab = ylab
  )
  graphics::axis(
    side = 2, at = grDevices::axisTicks(ylim, log = FALSE),
    las = ylas, cex.axis = ycex.axis
  )
  plotrix::plotCI(
    x = xseq, y = x, li = xlower, ui = xupper,
    pch = pch, cex = cex,
    sfrac = sfrac,
    col = col, add = TRUE
  )
  graphics::axis(side = 1, at = xseq, labels = NA)
  for (k in 1:length(xseq)) {
    graphics::axis(
      side = 1, at = xseq[k], labels = perf$names[k],
      col.axis = col.axis[k],
      las = xlas, cex.axis = xcex.axis
    )
  }
}

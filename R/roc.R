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
#' ids_train <- Resample(data = simul$ydata, tau = 0.5)
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' ids_recalib <- Resample(data = simul$ydata[-ids_train, , drop = FALSE], tau = 0.5)
#' xrecalib <- simul$xdata[ids_recalib, , drop = FALSE]
#' yrecalib <- simul$ydata[ids_recalib, , drop = FALSE]
#' xtest <- simul$xdata[-ids_recalib, ]
#' ytest <- simul$ydata[-ids_recalib, ]
#'
#' # Stability selection and recalibration
#' stab <- VariableSelection(xdata = xrecalib, ydata = yrecalib, family = "binomial")
#' recalibrated <- Recalibrate(xdata = xtrain, ydata = ytrain, stability = stab)
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
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' print(SelectedVariables(stab))
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   stability = stab
#' )
#' recalibrated$coefficients # recalibrated coefficients
#' head(recalibrated$fitted.values) # recalibrated predicted probabilities
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
#'
#' # Stability selection
#' ydata <- cbind(
#'   time = runif(nrow(simul$ydata), min = 100, max = 2000),
#'   case = simul$ydata[, 1]
#' ) # including dummy time to event
#' stab <- VariableSelection(xdata = simul$xdata, ydata = ydata, family = "cox")
#' print(SelectedVariables(stab))
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = simul$xdata, ydata = ydata,
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
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "binomial")
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   stability = stab
#' )
#' recalibrated$coefficients # recalibrated coefficients
#' head(recalibrated$fitted.values) # recalibrated predicted probabilities
#'
#'
#' ## Multinomial regression
#'
#' # Stability selection
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 15, family = "multinomial")
#' stab <- VariableSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "multinomial"
#' )
#'
#' # Recalibrating the model
#' recalibrated <- Recalibrate(
#'   xdata = simul$xdata, ydata = simul$ydata,
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
  myformula <- stats::as.formula(paste0("ydata ~ ", paste(ids, collapse = " + ")))

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


#' Prediction performance for regression models
#'
#' Computes the prediction performance for logistic or Cox regression models.
#' This is done by (i) recalibrating the model on a training set including a
#' proportion \code{tau} of the observations, and (ii) evaluating the
#' performance on the remaining observations (test set). For more reliable
#' results, the procedure is repeated \code{K} times (default \code{K=1}).
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
#' @family prediction performance functions
#'
#' @examples
#' \dontrun{
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 10, family = "binomial")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(data = simul$ydata, tau = 0.5)
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
#' ids_train <- Resample(data = simul$ydata, tau = 0.5)
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
                                   K = 1, tau = 0.7, seed = 1,
                                   n_thr = NULL,
                                   ij_method = FALSE, time = 3653) {
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
        ids_test <- Resample(data = ydata, tau = 1 - tau)
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


#' Receiver Operating Characteristic (ROC) curve
#'
#' Plots the True Positive Rate (TPR) as a function of the False Positive Rate
#' (FPR) for different thresholds in predicted probabilities. If the results
#' from multiple ROC analyses are provided (e.g. output from
#' \code{\link{ExplanatoryPerformance}} with large \code{K}), the point-wise
#' average is represented and flanked by a transparent band defined by
#' point-wise \code{quantiles}.
#'
#' @param roc output from \code{\link{ROC}} or
#'   \code{\link{ExplanatoryPerformance}}.
#' @param col colour of the point-wise average curve.
#' @param col_band colour of the band defined by point-wise \code{quantiles}.
#' @param alpha level of opacity for the band.
#' @param lwd line width for the point-wise average curve, as in
#'   \code{\link{par}}.
#' @param lty line type for the point-wise average curve, as in
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
#' ids_train <- Resample(data = simul$ydata, tau = 0.5)
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
PlotROC <- function(roc,
                    col = "red", col_band = "red",
                    alpha = 0.5,
                    lwd = 1, lty = 1,
                    quantiles = c(0.05, 0.95),
                    add = FALSE) {
  # Extracting the number of iterations
  niter <- length(roc$AUC)

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
  graphics::lines(apply(roc$FPR, 2, mean), apply(roc$TPR, 2, mean),
    type = "l", lwd = lwd, lty = lty, col = col
  )
}

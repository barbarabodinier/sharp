#' Regression model refitting
#'
#' Refits the regression model with stably selected variables as predictors
#' (without penalisation). Variables in \code{xdata} not evaluated in the
#' stability selection model will automatically be included as predictors.
#'
#' @inheritParams VariableSelection
#' @param stability output of \code{\link{VariableSelection}} or
#'   \code{\link{BiSelection}}. If \code{stability=NULL} (the default), a model
#'   including all variables in \code{xdata} as predictors is fitted. Argument
#'   \code{family} must be provided in this case.
#' @param family type of regression model. Possible values include
#'   \code{"gaussian"} (linear regression), \code{"binomial"} (logistic
#'   regression), \code{"multinomial"} (multinomial regression), and
#'   \code{"cox"} (survival analysis). If provided, this argument must be
#'   consistent with input \code{stability}.
#' @param implementation optional function to refit the model. If
#'   \code{implementation=NULL} and \code{stability} is the output of
#'   \code{\link{VariableSelection}}, \code{\link[stats]{lm}} (linear
#'   regression), \code{\link[survival]{coxph}} (Cox regression),
#'   \code{\link[stats]{glm}} (logistic regression), or
#'   \code{\link[nnet]{multinom}} (multinomial regression) is used. The function
#'   \code{\link{PLS}} is used for the output of \code{\link{BiSelection}}.
#' @param ... additional arguments to be passed to the function provided in
#'   \code{implementation}.
#'
#'
#' @return The output as obtained from: \item{\code{\link[stats]{lm}}}{for
#'   linear regression (\code{"gaussian"} family).}
#'   \item{\code{\link[survival]{coxph}}}{for Cox regression (\code{"cox"}
#'   family).} \item{\code{\link[stats]{glm}}}{for logistic regression
#'   (\code{"binomial"} family).} \item{\code{\link[nnet]{multinom}}}{for
#'   multinomial regression (\code{"multinomial"} family).}
#'
#' @seealso \code{\link{VariableSelection}}
#'
#' @examples
#' \donttest{
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
#' xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#' yrefit <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "gaussian")
#' print(SelectedVariables(stab))
#'
#' # Refitting the model
#' refitted <- Refit(
#'   xdata = xrefit, ydata = yrefit,
#'   stability = stab
#' )
#' refitted$coefficients # refitted coefficients
#' head(refitted$fitted.values) # refitted predicted values
#'
#' # Fitting the full model (including all possible predictors)
#' refitted <- Refit(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian"
#' )
#' refitted$coefficients # refitted coefficients
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
#' xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#' yrefit <- ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "cox")
#' print(SelectedVariables(stab))
#'
#' # Refitting the model
#' refitted <- Refit(
#'   xdata = xrefit, ydata = yrefit,
#'   stability = stab
#' )
#' refitted$coefficients # refitted coefficients
#' head(refitted$linear.predictors) # refitted scores
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
#' xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#' yrefit <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "binomial")
#'
#' # Refitting the model
#' refitted <- Refit(
#'   xdata = xrefit, ydata = yrefit,
#'   stability = stab
#' )
#' refitted$coefficients # refitted coefficients
#' head(refitted$fitted.values) # refitted predicted probabilities
#'
#'
#' ## Partial Least Squares (single component)
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
#' xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#' yrefit <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   implementation = SparsePLS,
#'   family = "gaussian"
#' )
#' print(SelectedVariables(stab))
#'
#' # Refitting the model
#' refitted <- Refit(
#'   xdata = xrefit, ydata = yrefit,
#'   implementation = PLS,
#'   stability = stab
#' )
#' refitted$Wmat # refitted X-weights
#' head(refitted$Tmat) # refitted X-scores
#'
#'
#' ## Partial Least Squares (multiple components)
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 15, q = 3, family = "gaussian")
#'
#' # Data split
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "gaussian"
#' )
#' xtrain <- simul$xdata[ids_train, , drop = FALSE]
#' ytrain <- simul$ydata[ids_train, , drop = FALSE]
#' xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#' yrefit <- simul$ydata[-ids_train, , drop = FALSE]
#'
#' # Stability selection
#' stab <- BiSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   family = "gaussian", ncomp = 3,
#'   LambdaX = 1:(ncol(xtrain) - 1),
#'   LambdaY = 1:(ncol(ytrain) - 1),
#'   implementation = SparsePLS
#' )
#' plot(stab)
#'
#' # Refitting the model
#' refitted <- Refit(
#'   xdata = xrefit, ydata = yrefit,
#'   stability = stab
#' )
#' refitted$Wmat # refitted X-weights
#' refitted$Cmat # refitted Y-weights
#' }
#'
#' @export
Refit <- function(xdata, ydata, stability = NULL,
                  family = NULL, implementation = NULL,
                  verbose = TRUE, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  if ("check_input" %in% names(extra_args)) {
    check_input <- extra_args$check_input
    extra_args <- extra_args[!names(extra_args) %in% "check_input"]
  } else {
    check_input <- TRUE
  }

  # Defining the type of model (PLS vs regression)
  use_pls <- FALSE

  # Checking input
  if (!is.null(stability)) {
    if (!inherits(stability, c("variable_selection", "bi_selection"))) {
      stop("Argument 'stability' is not of class 'variable_selection' or 'bi_selection'. This function can only be applied on the output of (i) VariableSelection() or (ii) BiSelection() for PLS models.")
    }
    if (inherits(stability, "bi_selection")) {
      # Checking mixOmics package is installed
      CheckPackageInstalled("mixOmics")
      use_pls <- TRUE

      if (!stability$methods$family %in% c("gaussian")) {
        stop("This function can only be applied with the 'gaussian' family for PLS models.")
      }
    } else {
      if (!stability$methods$family %in% c("gaussian", "cox", "binomial", "multinomial")) {
        stop("This function can only be applied with the following families for regression models: 'gaussian', 'cox', 'binomial', or 'multinomial'.")
      }
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

  # Object preparation, error and warning messages
  if (check_input) {
    CheckDataRegression(
      xdata = xdata, ydata = ydata, family = family, verbose = verbose
    )
  }

  # # Re-formatting the inputs
  # if (is.vector(xdata)) {
  #   xdata <- cbind(xdata)
  #   colnames(xdata) <- "var"
  # }
  # if (family %in% c("binomial", "multinomial")) {
  #   if (!is.factor(ydata)) {
  #     if (!is.vector(ydata)) {
  #       if (ncol(ydata) != 1) {
  #         ydata <- DummyToCategories(ydata)
  #       } else {
  #         ydata <- as.numeric(ydata)
  #       }
  #     }
  #   }
  # }

  if (use_pls) {
    # Refitting the PLS model
    mymodel <- PLS(
      xdata = xdata, ydata = ydata,
      selectedX = stability$selectedX,
      selectedY = stability$selectedY,
      family = family, ncomp = NULL,
      scale = stability$methods$scale
    )
  } else {
    # Extracting the stably selected predictors
    if (is.null(stability)) {
      selected <- rep(1, ncol(xdata))
      names(selected) <- colnames(xdata)
    } else {
      selected <- SelectedVariables(stability)
    }

    # Defining predictors for the model (including un-penalised)
    ids <- c(
      names(selected)[which(selected == 1)],
      colnames(xdata)[!colnames(xdata) %in% names(selected)]
    )

    if (is.null(implementation)) {
      # Writing model formula
      ids <- gsub("`", "", ids)
      colnames(xdata) <- gsub("`", "", colnames(xdata))
      if (length(ids) == 0) {
        message("No stably selected variables. Running a model with intercept only.")
        myformula <- stats::as.formula("ydata ~ 1")
      } else {
        myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", ids, "`"), collapse = " + ")))
      }

      # Recalibration for linear regression
      if (family == "gaussian") {
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::lm)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data")]
        mymodel <- do.call(stats::lm, args = c(
          list(
            formula = myformula,
            data = as.data.frame(xdata)
          ),
          tmp_extra_args
        ))
        # mymodel <- stats::lm(myformula, data = as.data.frame(xdata), ...)
      }

      # Recalibration for Cox regression
      if (family == "cox") {
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = survival::coxph)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data")]
        ydata <- survival::Surv(time = ydata[, 1], event = ydata[, 2])
        mymodel <- do.call(survival::coxph, args = c(
          list(
            formula = myformula,
            data = as.data.frame(xdata)
          ),
          tmp_extra_args
        ))
        # mymodel <- survival::coxph(myformula, data = as.data.frame(xdata), ...)
      }

      # Recalibration for logistic regression
      if (family == "binomial") {
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = stats::glm)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data", "family")]
        mymodel <- do.call(stats::glm, args = c(
          list(
            formula = myformula,
            data = as.data.frame(xdata),
            family = stats::binomial(link = "logit")
          ),
          tmp_extra_args
        ))
        # mymodel <- stats::glm(myformula,
        #   data = as.data.frame(xdata),
        #   family = stats::binomial(link = "logit"),
        #   ...
        # )
      }

      # Recalibration for multinomial regression
      if (family == "multinomial") {
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = nnet::multinom)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data", "trace")]
        mymodel <- do.call(nnet::multinom, args = c(
          list(
            formula = myformula,
            data = as.data.frame(xdata),
            trace = FALSE
          ),
          tmp_extra_args
        ))
        # mymodel <- nnet::multinom(myformula, data = as.data.frame(xdata), trace = FALSE, ...)
      }
    } else {
      tmp_extra_args <- extra_args[!names(extra_args) %in% c("xdata", "ydata", "family")]
      xdata <- xdata[, ids, drop = FALSE]
      mymodel <- do.call(implementation, args = c(
        list(
          xdata = xdata,
          ydata = ydata,
          family = family
        ),
        tmp_extra_args
      ))
    }
  }

  return(mymodel)
}


#' @rdname Refit
#' @export
Recalibrate <- Refit


#' Prediction performance in regression
#'
#' Calculates model performance for linear (measured by Q-squared), logistic
#' (AUC) or Cox (C-statistic) regression. This is done by (i) refitting the
#' model on a training set including a proportion \code{tau} of the
#' observations, and (ii) evaluating the performance on the remaining
#' observations (test set). For more reliable results, the procedure can be
#' repeated \code{K} times (default \code{K=1}).
#'
#' @inheritParams Refit
#' @param stability output of \code{\link{VariableSelection}}. If
#'   \code{stability=NULL} (the default), a model including all variables in
#'   \code{xdata} as predictors is fitted. Argument \code{family} must be
#'   provided in this case.
#' @param implementation optional function to refit the model. If
#'   \code{implementation=NULL} and \code{stability} is the output of
#'   \code{\link{VariableSelection}}, \code{\link[stats]{lm}} (linear
#'   regression), \code{\link[survival]{coxph}} (Cox regression),
#'   \code{\link[stats]{glm}} (logistic regression), or
#'   \code{\link[nnet]{multinom}} (multinomial regression) is used.
#' @param prediction optional function to compute predicted values from the
#'   model refitted with \code{implementation}.
#' @param K number of training-test splits.
#' @param tau proportion of observations used in the training set.
#' @param seed value of the seed to ensure reproducibility of the results.
#' @param n_thr number of thresholds to use to construct the ROC curve. If
#'   \code{n_thr=NULL}, all predicted probability values are iteratively used as
#'   thresholds. For faster computations on large data, less thresholds can be
#'   used. Only applicable to logistic regression.
#' @param time numeric indicating the time for which the survival probabilities
#'   are computed. Only applicable to Cox regression.
#' @param ij_method logical indicating if the analysis should be done for only
#'   one refitting/test split with variance of the concordance index should be
#'   computed using the infinitesimal jackknife method as implemented in
#'   \code{\link[survival]{concordance}}. If \code{ij_method=FALSE} (the
#'   default), the concordance indices computed for different refitting/test
#'   splits are reported. If \code{ij_method=TRUE}, the concordance index and
#'   estimated confidence interval at level 0.05 are reported. Only applicable
#'   to Cox regression.
#' @param resampling resampling approach to create the training set. The default
#'   is \code{"subsampling"} for sampling without replacement of a proportion
#'   \code{tau} of the observations. Alternatively, this argument can be a
#'   function to use for resampling. This function must use arguments named
#'   \code{data} and \code{tau} and return the IDs of observations to be
#'   included in the resampled dataset.
#' @param ... additional parameters passed to the function provided in
#'   \code{resampling}.
#'
#' @details For a fair evaluation of the prediction performance, the data is
#'   split into a training set (including a proportion \code{tau} of the
#'   observations) and test set (remaining observations). The regression model
#'   is fitted on the training set and applied on the test set. Performance
#'   metrics are computed in the test set by comparing predicted and observed
#'   outcomes.
#'
#'   For logistic regression, a Receiver Operating Characteristic (ROC) analysis
#'   is performed: the True and False Positive Rates (TPR and FPR), and Area
#'   Under the Curve (AUC) are computed for different thresholds in predicted
#'   probabilities.
#'
#'   For Cox regression, the Concordance Index (as implemented in
#'   \code{\link[survival]{concordance}}) looking at survival probabilities up
#'   to a specific \code{time} is computed.
#'
#'   For linear regression, the squared correlation between predicted and
#'   observed outcome in the test set (Q-squared) is reported.
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
#'   \code{ij_method=TRUE}).} \item{Beta}{matrix of estimated beta coefficients
#'   across the \code{K} iterations. Coefficients are extracted using the
#'   \code{\link[stats]{coef}} function.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Refit}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \donttest{
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 1000, pk = 10,
#'   family = "binomial", ev_xy = 0.7
#' )
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
#' # Evaluation of the performances on refitted models (K=1)
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, n_thr = NULL
#' )
#' plot(roc)
#'
#' # Using more refitting/test splits
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 100
#' )
#' boxplot(roc$AUC, ylab = "AUC")
#' plot(roc)
#'
#' # Comparison with saturated model
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   family = "binomial", K = 100
#' )
#' plot(roc, col = "blue", col_band = "blue", add = TRUE)
#'
#'
#' ## Partial Least Squares (single component)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   implementation = SparsePLS,
#'   family = "binomial"
#' )
#' print(SelectedVariables(stab))
#'
#' # Defining wrapping functions for PLS-DA
#' PLSDA <- function(xdata, ydata, family = "binomial") {
#'   model <- mixOmics::plsda(X = xdata, Y = as.factor(ydata), ncomp = 1)
#'   return(model)
#' }
#' PredictPLSDA <- function(xdata, model) {
#'   xdata <- xdata[, rownames(model$loadings$X), drop = FALSE]
#'   predicted <- predict(object = model, newdata = xdata)$predict[, 2, 1]
#'   return(predicted)
#' }
#'
#' # Evaluation of the performances on refitted models (K=1)
#' roc <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab,
#'   implementation = PLSDA, prediction = PredictPLSDA
#' )
#' plot(roc)
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
#' # Evaluation of the performances on refitted models (K=1)
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, ij_method = TRUE
#' )
#' print(perf)
#'
#' # Using more refitting/test splits
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab, K = 10, time = 1000
#' )
#' boxplot(perf$concordance)
#'
#'
#' ## Linear regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 10, family = "gaussian")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "gaussian"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "gaussian")
#'
#' # Evaluation of the performances on refitted models (K=1)
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab
#' )
#' print(perf)
#'
#'
#' ## Partial Least Squares (single component)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   implementation = SparsePLS,
#'   family = "gaussian"
#' )
#' print(SelectedVariables(stab))
#'
#' # Evaluation of the performances on refitted models (K=1)
#' perf <- ExplanatoryPerformance(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab,
#'   implementation = PLS, prediction = PredictPLS
#' )
#' print(perf)
#' }
#'
#' @export
ExplanatoryPerformance <- function(xdata, ydata,
                                   stability = NULL, family = NULL,
                                   implementation = NULL, prediction = NULL, resampling = "subsampling",
                                   K = 1, tau = 0.8, seed = 1,
                                   n_thr = NULL,
                                   ij_method = FALSE, time = 1000,
                                   verbose = TRUE, ...) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (!inherits(stability, "variable_selection")) {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!stability$methods$family %in% c("cox", "binomial", "gaussian")) {
      stop("This function can only be applied with the following families: 'binomial', 'cox' or 'gaussian'.")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'gaussian', 'cox' or 'binomial'.")
    }
  }

  # Object preparation, error and warning messages
  CheckDataRegression(
    xdata = xdata, ydata = ydata, family = family, verbose = verbose
  )

  # # Re-formatting input data
  # if (is.vector(ydata)) {
  #   ydata <- cbind(ydata)
  # }

  # Defining the metric to use
  if (ij_method) {
    K <- 1
  }
  if (family == "binomial") {
    metric <- "roc"
  }
  if (family == "cox") {
    metric <- "concordance"
  }
  if (family == "gaussian") {
    metric <- "q2"
  }
  # if (is.null(n_folds) | (K==1)) {
  #   n_folds <- 1
  # } else {
  #   K <- ceiling(K / n_folds)
  # }
  n_folds <- 1

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Defining the number of thresholds to use for AUC calculations
  if (is.null(n_thr)) {
    n_thr <- floor((1 - tau) * nrow(xdata))
  }

  # Running the subsampling iterations
  iter <- 0
  for (k in 1:K) {
    for (fold_id in 1:n_folds) {
      iter <- iter + 1
      if (n_folds == 1) {
        # Balanced training/test split
        ids_test <- Resample(data = ydata, tau = 1 - tau, family = family, resampling = resampling, ...)
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
      refitted <- Refit(
        xdata = xtrain, ydata = ytrain,
        stability = stability,
        implementation = implementation,
        family = family,
        check_input = FALSE,
        ...
      )

      if (is.null(implementation)) {
        # Initialising matrix of beta coefficients
        if (iter == 1) {
          if (family %in% c("gaussian", "binomial", "cox")) {
            Beta <- matrix(NA, nrow = K, ncol = length(stats::coef(refitted)))
            colnames(Beta) <- names(stats::coef(refitted))
            rownames(Beta) <- paste0("iter", 1:K)
          }
        }

        # Storing beta coefficients
        if (family %in% c("gaussian", "binomial", "cox")) {
          Beta[iter, ] <- stats::coef(refitted)
        }

        # Predictions from logistic models
        if (tolower(metric) == "roc") {
          predicted <- stats::predict.glm(refitted, newdata = as.data.frame(xtest), type = "response")
        }


        # Predictions from linear models
        if (tolower(metric) == "q2") {
          predicted <- stats::predict.lm(refitted, newdata = as.data.frame(xtest))
        }
      } else {
        if (is.null(prediction)) {
          stop("Argument 'prediction' has to be provided if 'implementation' is provided. It must be a function that takes the output of 'implementation' as argument.")
        }
        predicted <- do.call(prediction, args = list(xdata = xtest, model = refitted))
      }
      # Performing ROC analyses
      if (tolower(metric) == "roc") {
        # ROC analysis
        roc <- ROC(predicted = predicted, observed = ytest, n_thr = n_thr)

        # Initialisation of the object
        if (iter == 1) {
          # n_thr <- length(roc$FPR) - 2
          FPR <- TPR <- matrix(NA, nrow = K * n_folds, ncol = length(roc$TPR))
          AUC <- rep(NA, K * n_folds)
        }

        # Storing the metrics
        FPR[iter, ] <- roc$FPR
        TPR[iter, ] <- roc$TPR
        AUC[iter] <- roc$AUC
      }

      # Performing Q-squared analyses
      if (tolower(metric) == "q2") {
        # Initialisation of the object
        if (iter == 1) {
          Q_squared <- rep(NA, K * n_folds)
        }

        # Computing the Q-squared
        Q_squared[iter] <- stats::cor(predicted, ytest)^2
      }

      # Performing concordance analyses
      if (tolower(metric) == "concordance") {
        # Computing the concordance index for given times
        predicted <- stats::predict(refitted, newdata = as.data.frame(xtest), type = "lp")
        survobject <- survival::Surv(time = ytest[, 1], event = ytest[, 2])
        S0 <- summary(survival::survfit(refitted), times = time, extend = TRUE)$surv
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

  # Preparing the output
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

  if (tolower(metric) == "q2") {
    out <- list(Q_squared = Q_squared)
  }

  if (is.null(implementation)) {
    out <- c(out, Beta = list(Beta))
  }

  # Defining class
  class(out) <- "roc_band"

  return(out)
}


#' Incremental prediction performance in regression
#'
#' Computes the prediction performance of regression models where predictors are
#' sequentially added by order of decreasing selection proportion. This function
#' can be used to evaluate the marginal contribution of each of the selected
#' predictors over and above more stable predictors. Performances are evaluated
#' as in \code{\link{ExplanatoryPerformance}}.
#'
#' @inheritParams ExplanatoryPerformance
#' @param n_predictors number of predictors to consider.
#' @param verbose logical indicating if a loading bar and messages should be
#'   printed.
#'
#' @return An object of class \code{incremental}.
#'
#'   For logistic regression, a list with: \item{FPR}{A list with, for each of
#'   the models (sequentially added predictors), the False Positive Rates for
#'   different thresholds (columns) and different data splits (rows).}
#'   \item{TPR}{A list with, for each of the models (sequentially added
#'   predictors), the True Positive Rates for different thresholds (columns) and
#'   different data splits (rows).} \item{AUC}{A list with, for each of the
#'   models (sequentially added predictors), a vector of Area Under the Curve
#'   (AUC) values obtained with different data splits.} \item{Beta}{Estimated
#'   regression coefficients from visited models.} \item{names}{Names of the
#'   predictors by order of inclusion.}
#'
#'   For Cox regression, a list with: \item{concordance}{If
#'   \code{ij_method=FALSE}, a list with, for each of the models (sequentially
#'   added predictors), a vector of concordance indices obtained with different
#'   data splits. If \code{ij_method=TRUE}, a vector of concordance indices for
#'   each of the models (sequentially added predictors).} \item{lower}{A vector
#'   of the lower bound of the confidence interval at level 0.05 for concordance
#'   indices for each of the models (sequentially added predictors). Only
#'   returned if \code{ij_method=TRUE}.} \item{upper}{A vector of the upper
#'   bound of the confidence interval at level 0.05 for concordance indices for
#'   each of the models (sequentially added predictors). Only returned if
#'   \code{ij_method=TRUE}.} \item{Beta}{Estimated regression coefficients from
#'   visited models.} \item{names}{Names of the predictors by order of
#'   inclusion.}
#'
#'   For linear regression, a list with: \item{Q_squared}{A list with, for each
#'   of the models (sequentially added predictors), a vector of Q-squared
#'   obtained with different data splits.} \item{Beta}{Estimated regression
#'   coefficients from visited models.} \item{names}{Names of the predictors by
#'   order of inclusion.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Refit}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \donttest{
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
#' summary(perf)
#'
#' # Visualisation
#' PlotIncremental(perf)
#' plot(perf) # alternative formulation
#'
#'
#' ## Partial Least Squares (single component)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   implementation = SparsePLS,
#'   family = "binomial"
#' )
#' print(SelectedVariables(stab))
#'
#' # Defining wrapping functions for PLS-DA
#' PLSDA <- function(xdata, ydata, family = "binomial") {
#'   model <- mixOmics::plsda(X = xdata, Y = as.factor(ydata), ncomp = 1)
#'   return(model)
#' }
#' PredictPLSDA <- function(xdata, model) {
#'   xdata <- xdata[, rownames(model$loadings$X), drop = FALSE]
#'   predicted <- predict(object = model, newdata = xdata)$predict[, 2, 1]
#'   return(predicted)
#' }
#'
#' # Evaluation of the performances on refitted models (K=1)
#' incremental <- Incremental(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab,
#'   implementation = PLSDA, prediction = PredictPLSDA,
#'   K = 10
#' )
#' PlotIncremental(incremental)
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
#'
#'
#' ## Linear regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "gaussian")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "gaussian"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "gaussian")
#'
#' # Evaluating marginal contribution of the predictors
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#' PlotIncremental(perf)
#'
#'
#' ## Partial Least Squares (single component)
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xtrain, ydata = ytrain,
#'   implementation = SparsePLS,
#'   family = "gaussian"
#' )
#' print(SelectedVariables(stab))
#'
#' # Evaluation of the performances on refitted models (K=1)
#' incremental <- Incremental(
#'   xdata = xtest, ydata = ytest,
#'   stability = stab,
#'   implementation = PLS, prediction = PredictPLS,
#'   K = 10
#' )
#' PlotIncremental(incremental)
#' }
#'
#' @export
Incremental <- function(xdata, ydata,
                        stability = NULL, family = NULL,
                        implementation = NULL, prediction = NULL, resampling = "subsampling",
                        n_predictors = NULL,
                        K = 100, tau = 0.8, seed = 1,
                        n_thr = NULL,
                        ij_method = FALSE, time = 1000,
                        verbose = TRUE, ...) {
  # Checking the inputs
  if (!is.null(stability)) {
    if (!inherits(stability, "variable_selection")) {
      stop("Argument 'stability' is not of class 'variable_selection'. This function can only be applied on the output of VariableSelection().")
    }
    if (!is.null(family)) {
      if (family != stability$methods$family) {
        warning(paste0("Arguments 'stability' and 'family' are not consistent. The family specified in argument stability was used: ", stability$methods$family))
      }
    }
    family <- stability$methods$family
    if (!family %in% c("cox", "binomial", "gaussian")) {
      stop("This function can only be applied with the following families: 'binomial', 'cox' or 'gaussian'.")
    }
  } else {
    if (is.null(family)) {
      stop("Argument 'family' must be provided. Possible values are: 'binomial', 'cox' or 'gaussian'.")
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
  Beta <- list()
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
  if (family == "gaussian") {
    Q_squared <- list()
  }

  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }

  for (k in 1:n_predictors) {
    perf <- ExplanatoryPerformance(
      xdata = xdata[, myorder[1:k], drop = FALSE],
      ydata = ydata,
      stability = NULL,
      family = family,
      implementation = implementation,
      prediction = prediction,
      resampling = resampling,
      K = K, tau = tau, seed = seed,
      n_thr = n_thr,
      ij_method = ij_method, time = time,
      ...
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
    if (family == "gaussian") {
      Q_squared <- c(Q_squared, list(perf$Q_squared))
    }
    Beta <- c(Beta, list(perf$Beta))

    if (verbose) {
      utils::setTxtProgressBar(pb, k / n_predictors)
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
  if (family == "gaussian") {
    out <- list(Q_squared = Q_squared)
  }

  # Adding beta coefficients
  if (is.null(implementation)) {
    out <- c(out, list(Beta = Beta))
  }

  # Adding variable names
  out <- c(out, names = list(myorder[1:n_predictors]))

  # Defining class
  class(out) <- "incremental"

  return(out)
}


#' Visualisation of incremental performance
#'
#' Represents prediction performances upon sequential inclusion of the
#' predictors in a logistic or Cox regression model as produced by
#' \code{\link{Incremental}}. The median and \code{quantiles} of the performance
#' metric are reported.
#'
#' @inheritParams CalibrationPlot
#' @param perf output of \code{\link{Incremental}}.
#' @param quantiles quantiles defining the lower and upper bounds.
#' @param sfrac size of the end bars, as in \code{\link[plotrix]{plotCI}}.
#' @param col vector of point colours.
#' @param col.axis optional vector of label colours. If \code{col.axis=NULL},
#'   the colours provided in argument \code{col} are used.
#' @param xcex.axis size of labels along the x-axis.
#' @param ycex.axis size of labels along the y-axis.
#' @param output_data logical indicating if the median and quantiles should be
#'   returned in a matrix.
#'
#' @return A plot.
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Refit}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \donttest{
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
#'
#'
#' ## Linear regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "gaussian")
#'
#' # Balanced split: 50% variable selection set and 50% for evaluation of performances
#' ids_train <- Resample(
#'   data = simul$ydata,
#'   tau = 0.5, family = "gaussian"
#' )
#' xtrain <- simul$xdata[ids_train, ]
#' ytrain <- simul$ydata[ids_train, ]
#' xtest <- simul$xdata[-ids_train, ]
#' ytest <- simul$ydata[-ids_train, ]
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = xtrain, ydata = ytrain, family = "gaussian")
#'
#' # Evaluating marginal contribution of the predictors
#' perf <- Incremental(xdata = xtest, ydata = ytest, stability = stab, K = 10)
#' PlotIncremental(perf)
#'
#' # Evaluating marginal contribution of the predictors beyond stably selected
#' perf <- Incremental(
#'   xdata = xtest, ydata = ytest, stability = stab,
#'   K = 10, n_predictors = 20
#' )
#' N <- sum(SelectedVariables(stab))
#' PlotIncremental(perf, col = c(rep("red", N), rep("grey", 20 - N)))
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
                            xgrid = FALSE, ygrid = FALSE,
                            output_data = FALSE) {
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
      x <- sapply(perf$concordance, stats::median, na.rm = TRUE)
      xlower <- sapply(perf$concordance, stats::quantile, probs = quantiles[1], na.rm = TRUE)
      xupper <- sapply(perf$concordance, stats::quantile, probs = quantiles[2], na.rm = TRUE)
    }
  }

  if ("AUC" %in% names(perf)) {
    x <- sapply(perf$AUC, stats::median, na.rm = TRUE)
    xlower <- sapply(perf$AUC, stats::quantile, probs = quantiles[1], na.rm = TRUE)
    xupper <- sapply(perf$AUC, stats::quantile, probs = quantiles[2], na.rm = TRUE)
  }

  if ("Q_squared" %in% names(perf)) {
    x <- sapply(perf$Q_squared, stats::median, na.rm = TRUE)
    xlower <- sapply(perf$Q_squared, stats::quantile, probs = quantiles[1], na.rm = TRUE)
    xupper <- sapply(perf$Q_squared, stats::quantile, probs = quantiles[2], na.rm = TRUE)
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
      side = 1, at = xseq[k],
      labels = ifelse(k == 1, yes = perf$names[k], no = paste0("+ ", perf$names[k])),
      col.axis = col.axis[k],
      las = xlas, cex.axis = xcex.axis
    )
  }

  if (output_data) {
    mat <- rbind(x, xlower, xupper)
    colnames(mat) <- perf$names
    return(mat)
  }
}

#' @rdname PlotIncremental
#' @export
IncrementalPlot <- PlotIncremental

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
#' @param implementation optional function to refit the model. If \code{stability} is the output of
#'   \code{\link{VariableSelection}}, a regression model is refitted.
#'   If \code{implementation=NULL} and \code{Lambda=0}, this is done using \code{\link[stats]{lm}} (for linear
#'   regression), \code{\link[survival]{coxph}} (Cox regression),
#'   \code{\link[stats]{glm}} (logistic regression), or
#'   \code{\link[nnet]{multinom}} (multinomial regression).
#'   If \code{Lambda=NULL}, a Ridge regression is fitted and calibrated by cross validation using \code{\link[glmnet]{cv.glmnet}}.
#'   The function \code{\link{PLS}} is used if \code{stability} is the output of \code{\link{BiSelection}}.
#' @param Lambda optional vector of penalty parameters.
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
#' ## Partial Least Squares (multiple components)
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 500, pk = 15, q = 3, family = "gaussian")
#'
#'   # Data split
#'   ids_train <- Resample(
#'     data = simul$ydata,
#'     tau = 0.5, family = "gaussian"
#'   )
#'   xtrain <- simul$xdata[ids_train, , drop = FALSE]
#'   ytrain <- simul$ydata[ids_train, , drop = FALSE]
#'   xrefit <- simul$xdata[-ids_train, , drop = FALSE]
#'   yrefit <- simul$ydata[-ids_train, , drop = FALSE]
#'
#'   # Stability selection
#'   stab <- BiSelection(
#'     xdata = xtrain, ydata = ytrain,
#'     family = "gaussian", ncomp = 3,
#'     LambdaX = seq_len(ncol(xtrain) - 1),
#'     LambdaY = seq_len(ncol(ytrain) - 1),
#'     implementation = SparsePLS
#'   )
#'   plot(stab)
#'
#'   # Refitting the model
#'   refitted <- Refit(
#'     xdata = xrefit, ydata = yrefit,
#'     stability = stab
#'   )
#'   refitted$Wmat # refitted X-weights
#'   refitted$Cmat # refitted Y-weights
#' }
#' }
#'
#' @export
Refit <- function(xdata, ydata, stability = NULL,
                  family = NULL, implementation = NULL,
                  Lambda = NULL, seed = 1,
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
    xdata <- xdata[, ids, drop = FALSE]

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

      # Defining penalisation
      penalised <- TRUE
      if (ncol(xdata) == 1) {
        penalised <- FALSE
      }

      # Preparing the model
      if (penalised) {
        # Extracting relevant extra arguments
        tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = glmnet::cv.glmnet)
        tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("x", "y", "lambda", "alpha", "family", "type.multinomial")]

        # Setting seed for reproducibility
        withr::local_seed(seed)

        # Running model
        if ((is.null(Lambda)) | (length(Lambda) > 1)) {
          mymodel <- do.call(glmnet::cv.glmnet, args = c(
            list(
              x = xdata,
              y = ydata,
              family = family,
              alpha = 0,
              lambda = Lambda,
              type.multinomial = "grouped"
            ),
            tmp_extra_args
          ))
        } else {
          mymodel <- do.call(glmnet::glmnet, args = c(
            list(
              x = xdata,
              y = ydata,
              family = family,
              alpha = 0,
              lambda = Lambda,
              type.multinomial = "grouped"
            ),
            tmp_extra_args
          ))
        }
      } else {
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
        }

        # Recalibration for Cox regression
        if (family == "cox") {
          ids1 <- which(names(extra_args) %in% names(formals(survival::coxph)))
          ids2 <- which(names(extra_args) %in% names(formals(survival::coxph.control)))
          tmp_extra_args <- extra_args[unique(c(ids1, ids2))]
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data")]
          ydata <- survival::Surv(time = ydata[, 1], event = ydata[, 2])
          mymodel <- do.call(survival::coxph, args = c(
            list(
              formula = myformula,
              data = as.data.frame(xdata)
            ),
            tmp_extra_args
          ))
        }

        # Recalibration for logistic regression
        if (family == "binomial") {
          ids1 <- which(names(extra_args) %in% names(formals(stats::glm)))
          ids2 <- which(names(extra_args) %in% names(formals(stats::glm.control)))
          tmp_extra_args <- extra_args[unique(c(ids1, ids2))]
          tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("formula", "data", "family")]
          suppressWarnings({
            mymodel <- do.call(stats::glm, args = c(
              list(
                formula = myformula,
                data = as.data.frame(xdata),
                family = stats::binomial(link = "logit")
              ),
              tmp_extra_args
            ))
          })
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
        }
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
#' @param family type of regression model. Possible values include
#'   \code{"gaussian"} (linear regression), \code{"binomial"} (logistic
#'   regression), and \code{"cox"} (survival analysis). If provided, this
#'   argument must be consistent with input \code{stability}.
#' @param new_xdata optional test set (predictor data).
#' @param new_ydata optional test set (outcome data).
#' @param implementation optional function to refit the model. If
#'   \code{implementation=NULL} and \code{stability} is the output of
#'   \code{\link{VariableSelection}}, \code{\link[stats]{lm}} (linear
#'   regression), \code{\link[survival]{coxph}} (Cox regression),
#'   \code{\link[stats]{glm}} (logistic regression), or
#'   \code{\link[nnet]{multinom}} (multinomial regression) is used.
#' @param prediction optional function to compute predicted values from the
#'   model refitted with \code{implementation}.
#' @param K number of training-test splits. Only used if \code{new_xdata} and
#'   \code{new_ydata} are not provided.
#' @param tau proportion of observations used in the training set. Only used if
#'   \code{new_xdata} and \code{new_ydata} are not provided.
#' @param seed value of the seed to ensure reproducibility of the results. Only
#'   used if \code{new_xdata} and \code{new_ydata} are not provided.
#' @param n_thr number of thresholds to use to construct the ROC curve. If
#'   \code{n_thr=NULL}, all predicted probability values are iteratively used as
#'   thresholds. For faster computations on large data, less thresholds can be
#'   used. Only applicable to logistic regression.
#' @param time numeric indicating the time for which the survival probabilities
#'   are computed. Only applicable to Cox regression.
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
#'   \item{Beta}{matrix of estimated beta coefficients across the \code{K}
#'   iterations. Coefficients are extracted using the \code{\link[stats]{coef}}
#'   function.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Refit}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 1000, pk = 20,
#'   family = "binomial", ev_xy = 0.8
#' )
#'
#' # Data split: selection, training and test set
#' ids <- Split(
#'   data = simul$ydata,
#'   family = "binomial",
#'   tau = c(0.4, 0.3, 0.3)
#' )
#' xselect <- simul$xdata[ids[[1]], ]
#' yselect <- simul$ydata[ids[[1]], ]
#' xtrain <- simul$xdata[ids[[2]], ]
#' ytrain <- simul$ydata[ids[[2]], ]
#' xtest <- simul$xdata[ids[[3]], ]
#' ytest <- simul$ydata[ids[[3]], ]
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xselect,
#'   ydata = yselect,
#'   family = "binomial"
#' )
#'
#' # Performances in test set of model refitted in training set
#' roc <- ExplanatoryPerformance(
#'   xdata = xtrain, ydata = ytrain,
#'   new_xdata = xtest, new_ydata = ytest,
#'   stability = stab
#' )
#' plot(roc)
#' roc$AUC
#'
#' # Alternative with multiple training/test splits
#' roc <- ExplanatoryPerformance(
#'   xdata = rbind(xtrain, xtest),
#'   ydata = c(ytrain, ytest),
#'   stability = stab, K = 100
#' )
#' plot(roc)
#' boxplot(roc$AUC)
#'
#' # Partial Least Squares Discriminant Analysis
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   stab <- VariableSelection(
#'     xdata = xselect,
#'     ydata = yselect,
#'     implementation = SparsePLS,
#'     family = "binomial"
#'   )
#'
#'   # Defining wrapping functions for predictions from PLS-DA
#'   PLSDA <- function(xdata, ydata, family = "binomial") {
#'     model <- mixOmics::plsda(X = xdata, Y = as.factor(ydata), ncomp = 1)
#'     return(model)
#'   }
#'   PredictPLSDA <- function(xdata, model) {
#'     xdata <- xdata[, rownames(model$loadings$X), drop = FALSE]
#'     predicted <- predict(object = model, newdata = xdata)$predict[, 2, 1]
#'     return(predicted)
#'   }
#'
#'   # Performances with custom models
#'   roc <- ExplanatoryPerformance(
#'     xdata = rbind(xtrain, xtest),
#'     ydata = c(ytrain, ytest),
#'     stability = stab, K = 100,
#'     implementation = PLSDA, prediction = PredictPLSDA
#'   )
#'   plot(roc)
#' }
#' }
#' @export
ExplanatoryPerformance <- function(xdata, ydata, new_xdata = NULL, new_ydata = NULL,
                                   stability = NULL, family = NULL,
                                   implementation = NULL, prediction = NULL, resampling = "subsampling",
                                   K = 1, tau = 0.8, seed = 1,
                                   n_thr = NULL,
                                   time = 1000,
                                   verbose = FALSE, ...) {
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
  if (!is.null(new_xdata)) {
    if (is.null(new_ydata)) {
      stop("Argument 'new_ydata' must be provided if 'new_xdata' is provided.")
    }
    K <- 1
  }

  # Object preparation, error and warning messages
  CheckDataRegression(
    xdata = xdata, ydata = ydata, family = family, verbose = verbose
  )
  if (!is.null(new_xdata)) {
    xtrain <- xdata
    ytrain <- ydata
    CheckDataRegression(
      xdata = new_xdata, ydata = new_ydata, family = family, verbose = verbose
    )
    xtest <- xdata # variable updated by CheckDataRegression()
    ytest <- ydata
  }

  # Defining the metric to use
  if (family == "binomial") {
    metric <- "roc"
  }
  if (family == "cox") {
    metric <- "concordance"
  }
  if (family == "gaussian") {
    metric <- "q2"
  }

  # Setting seed for reproducibility
  withr::local_seed(seed)

  # Running the subsampling iterations
  iter <- 0
  for (k in seq_len(K)) {
    iter <- iter + 1

    if (is.null(new_xdata)) {
      # Balanced training/test split
      ids_test <- Resample(data = ydata, tau = 1 - tau, family = family, resampling = resampling, ...)
      xtrain <- xdata[-ids_test, , drop = FALSE]
      ytrain <- ydata[-ids_test, , drop = FALSE]
      xtest <- xdata[ids_test, , drop = FALSE]
      ytest <- ydata[ids_test, , drop = FALSE]
    }

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
          tmpbeta <- as.vector(stats::coef(refitted))
          Beta <- matrix(NA, nrow = K, ncol = length(tmpbeta))
          colnames(Beta) <- names(tmpbeta)
          rownames(Beta) <- paste0("iter", seq_len(K))
        }
      }

      # Storing beta coefficients
      if (family %in% c("gaussian", "binomial", "cox")) {
        Beta[iter, ] <- as.vector(stats::coef(refitted))
      }

      if (inherits(refitted, c("glm", "lm"))) {
        # if (ncol(xtrain) == 1) {
        # Predictions from logistic models
        if (tolower(metric) == "roc") {
          suppressWarnings({
            predicted <- stats::predict.glm(refitted, newdata = as.data.frame(xtest))
          })
        }

        # Predictions from linear models
        if (tolower(metric) == "q2") {
          predicted <- stats::predict.lm(refitted, newdata = as.data.frame(xtest))
        }
      } else {
        ids_predictors <- rownames(stats::coef(refitted))
        ids_predictors <- ids_predictors[which(ids_predictors %in% colnames(xtest))]
        predicted <- stats::predict(refitted, newx = as.matrix(xtest[, ids_predictors, drop = FALSE]))
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
      roc <- ROC(predicted = predicted, observed = ytest, n_thr = length(predicted) - 1)

      # Initialisation of the object
      if (iter == 1) {
        FPR <- TPR <- matrix(NA, nrow = K, ncol = length(roc$TPR))
        AUC <- rep(NA, K)
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
        Q_squared <- rep(NA, K)
      }

      # Computing the Q-squared
      Q_squared[iter] <- stats::cor(predicted, ytest)^2
    }

    # Performing concordance analyses
    if (tolower(metric) == "concordance") {
      # Computing the concordance index for given times
      if (length(stats::coef(refitted)) == 1) {
        predicted <- stats::predict(refitted, newdata = as.data.frame(xtest), type = "lp")
      } else {
        predicted <- stats::predict(refitted, newx = as.matrix(xtest), type = "link")
      }
      ytest[which(ytest[, "time"] > time), "status"] <- 0
      cstat <- glmnet::Cindex(pred = predicted, y = ytest)
      if (iter == 1) {
        cindex <- rep(NA, K)
      }
      cindex[iter] <- cstat
    }
  }

  # Preparing the output
  if (tolower(metric) == "roc") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }

  if (tolower(metric) == "concordance") {
    out <- list(concordance = cindex)
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
#'   predictors by order of inclusion.} \item{stable}{Binary vector indicating
#'   which predictors are stably selected. Only returned if \code{stability} is
#'   provided.}
#'
#'   For Cox regression, a list with: \item{concordance}{A list with, for each
#'   of the models (sequentially added predictors), a vector of concordance
#'   indices obtained with different data splits.} \item{Beta}{Estimated
#'   regression coefficients from visited models.} \item{names}{Names of the
#'   predictors by order of inclusion.} \item{stable}{Binary vector indicating
#'   which predictors are stably selected. Only returned if \code{stability} is
#'   provided.}
#'
#'   For linear regression, a list with: \item{Q_squared}{A list with, for each
#'   of the models (sequentially added predictors), a vector of Q-squared
#'   obtained with different data splits.} \item{Beta}{Estimated regression
#'   coefficients from visited models.} \item{names}{Names of the predictors by
#'   order of inclusion.} \item{stable}{Binary vector indicating which
#'   predictors are stably selected. Only returned if \code{stability} is
#'   provided.}
#'
#' @seealso \code{\link{VariableSelection}}, \code{\link{Refit}}
#'
#' @family prediction performance functions
#'
#' @examples
#' \donttest{
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(
#'   n = 1000, pk = 20,
#'   family = "binomial", ev_xy = 0.8
#' )
#'
#' # Data split: selection, training and test set
#' ids <- Split(
#'   data = simul$ydata,
#'   family = "binomial",
#'   tau = c(0.4, 0.3, 0.3)
#' )
#' xselect <- simul$xdata[ids[[1]], ]
#' yselect <- simul$ydata[ids[[1]], ]
#' xtrain <- simul$xdata[ids[[2]], ]
#' ytrain <- simul$ydata[ids[[2]], ]
#' xtest <- simul$xdata[ids[[3]], ]
#' ytest <- simul$ydata[ids[[3]], ]
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = xselect,
#'   ydata = yselect,
#'   family = "binomial"
#' )
#'
#' # Performances in test set of model refitted in training set
#' incr <- Incremental(
#'   xdata = xtrain, ydata = ytrain,
#'   new_xdata = xtest, new_ydata = ytest,
#'   stability = stab, n_predictors = 10
#' )
#' plot(incr)
#'
#' # Alternative with multiple training/test splits
#' incr <- Incremental(
#'   xdata = rbind(xtrain, xtest),
#'   ydata = c(ytrain, ytest),
#'   stability = stab, K = 10, n_predictors = 10
#' )
#' plot(incr)
#' }
#'
#' @export
Incremental <- function(xdata, ydata, new_xdata = NULL, new_ydata = NULL,
                        stability = NULL, family = NULL,
                        implementation = NULL, prediction = NULL, resampling = "subsampling",
                        n_predictors = NULL,
                        K = 100, tau = 0.8, seed = 1,
                        n_thr = NULL,
                        time = 1000,
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

  # Object preparation, error and warning messages
  CheckDataRegression(
    xdata = xdata, ydata = ydata, family = family, verbose = verbose
  )

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

  # Defining the order of inclusion in the model
  if (is.null(stability)) {
    myorder <- colnames(xdata) # order of the columns is used
    n_predictors <- min(n_predictors, ncol(xdata))
  } else {
    # Identifying unpenalised variables
    unpenalised_vars <- colnames(xdata)[!colnames(xdata) %in% names(SelectedVariables(stability))]
    n_unpenalised <- length(unpenalised_vars)
    n_predictors <- min(n_predictors, ncol(xdata) - n_unpenalised)

    # Including variables by order of decreasing selection proportions
    myorder <- names(SelectionProportions(stability))[sort.list(SelectionProportions(stability), decreasing = TRUE)]
    myorder <- myorder[seq_len(n_predictors)]

    # Including the variables that are forced in the model first by order of columns in the data
    myorder <- c(unpenalised_vars, myorder)
  }

  # Initialisation of the objects
  Beta <- list()
  if (family == "binomial") {
    TPR <- FPR <- AUC <- list()
  }
  if (family == "cox") {
    concordance <- list()
  }
  if (family == "gaussian") {
    Q_squared <- list()
  }

  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }

  for (k in seq_len(length(myorder))) {
    perf <- ExplanatoryPerformance(
      xdata = xdata[, myorder[seq_len(k)], drop = FALSE],
      ydata = ydata,
      new_xdata = new_xdata,
      new_ydata = new_ydata,
      stability = NULL,
      family = family,
      implementation = implementation,
      prediction = prediction,
      resampling = resampling,
      K = K, tau = tau, seed = seed,
      n_thr = n_thr,
      time = time,
      verbose = FALSE,
      ...
    )
    if (family == "binomial") {
      FPR <- c(FPR, list(perf$FPR))
      TPR <- c(TPR, list(perf$TPR))
      AUC <- c(AUC, list(perf$AUC))
    }
    if (family == "cox") {
      concordance <- c(concordance, list(perf$concordance))
    }
    if (family == "gaussian") {
      Q_squared <- c(Q_squared, list(perf$Q_squared))
    }
    Beta <- c(Beta, list(perf$Beta))

    if (verbose) {
      utils::setTxtProgressBar(pb, k / length(myorder))
    }
  }

  # Preparing the output
  if (family == "binomial") {
    out <- list(FPR = FPR, TPR = TPR, AUC = AUC)
  }
  if (family == "cox") {
    out <- list(concordance = concordance)
  }
  if (family == "gaussian") {
    out <- list(Q_squared = Q_squared)
  }

  # Adding beta coefficients
  if (is.null(implementation)) {
    out <- c(out, list(Beta = Beta))
  }

  # Adding variable names
  mynames <- myorder
  out <- c(out, names = list(mynames))

  # Adding stable selection status
  if (!is.null(stability)) {
    mystable <- rep(0, length(myorder))
    names(mystable) <- myorder
    mystable[seq_len(max(which(myorder %in% names(which(SelectedVariables(stability) == 1)))))] <- 1
    out <- c(out, stable = list(mystable))
  }

  # Defining class
  class(out) <- "incremental"

  return(out)
}

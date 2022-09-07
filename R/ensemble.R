#' Ensemble model
#'
#' Creates an ensemble predictive model from \code{\link{VariableSelection}}
#' outputs.
#'
#' @inheritParams Refit
#' @param stability output of \code{\link{VariableSelection}}.
#'
#' @return An object of class \code{ensemble_model}. A list with:
#'   \item{intercept}{a vector of refitted intercepts for the \code{K}
#'   calibrated models.} \item{beta}{a matrix of beta coefficients from the
#'   \code{K} calibrated models.} \item{models}{a list of \code{K} models that
#'   can be used for prediction. These models are of class \code{"lm"} if
#'   \code{family="gaussian"} or \code{"glm"} if
#'   \code{family="binomial"}.} \item{family}{type of regression, extracted from
#'   \code{stability}. Possible values are \code{"gaussian"} or
#'   \code{"binomial"}.}
#'
#' @family ensemble model functions
#'
#' @examples
#' # Linear regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "gaussian")
#' ensemble <- Ensemble(stability = stab, xdata = simul$xdata, ydata = simul$ydata)
#'
#' # Logistic regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 200, pk = 20, family = "binomial")
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, family = "binomial")
#' ensemble <- Ensemble(stability = stab, xdata = simul$xdata, ydata = simul$ydata)
#'
#' @export
Ensemble <- function(stability, xdata, ydata) {
  # Checking family argument
  if (!stability$methods$family %in% c("gaussian", "binomial")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian' or 'binomial'.")
  }

  # Extracting beta coefficients
  argmax_id <- ArgmaxId(stability = stability)[1]
  beta <- t(stability$Beta[argmax_id, , ])

  # Checking xdata input
  if (!all(colnames(xdata) %in% colnames(beta))) {
    stop("Arguments 'stability' and 'xdata' are not consistent. Column names in 'xdata' do not match the ones used to create 'stability'.")
  }
  xdata <- xdata[, colnames(beta)]

  # Defining regression formula
  myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", colnames(xdata), "`"), collapse = " + ")))

  # Linear regression
  if (stability$methods$family == "gaussian") {
    # Checking ydata input
    if (inherits(ydata, "matrix")) {
      if (ncol(ydata) > 1) {
        stop("Arguments 'stability' and 'ydata' are not consistent. Matrix 'ydata' should only have one column.")
      }
    }

    # Creating lm template
    template <- stats::lm(myformula, data = as.data.frame(xdata))

    # Initialising objects
    intercept <- rep(NA, nrow(beta))
    models <- list()
    for (k in 1:nrow(beta)) {
      # Calculating intercept for specific beta coefficients
      b <- cbind(beta[k, ])
      # a=mean(ydata)-apply(xdata, 2, mean)%*%b
      tmpy <- ydata - as.matrix(xdata) %*% b
      a <- stats::coef(stats::lm(tmpy ~ 1))
      intercept[k] <- a

      # Setting coefficients in lm model
      template$coefficients <- c(`(Intercept)` = a, b)
      models <- c(models, list(template))
    }
  }

  # Logistic regression
  if (stability$methods$family == "binomial") {
    # Checking ydata input
    if (inherits(ydata, "matrix")) {
      if (ncol(ydata) > 1) {
        stop("Arguments 'stability' and 'ydata' are not consistent. Matrix 'ydata' should only have one column.")
      }
    }

    # Creating lm template
    template <- stats::glm(myformula,
      data = as.data.frame(xdata),
      family = stats::binomial(link = "logit")
    )

    # Initialising objects
    intercept <- rep(NA, nrow(beta))
    models <- list()
    for (k in 1:nrow(beta)) {
      # Calculating intercept for specific beta coefficients
      b <- cbind(beta[k, ])
      tmp <- as.vector(as.matrix(xdata) %*% b)
      a <- unname(stats::coef(stats::glm(ydata ~ 1, offset = tmp, family = stats::binomial(link = "logit"))))
      intercept[k] <- a

      # Setting coefficients in lm model
      template$coefficients <- c(`(Intercept)` = a, b)
      models <- c(models, list(template))
    }
  }

  # Preparing output
  out <- list(
    intercept = intercept,
    beta = beta,
    models = models,
    family = stability$methods$family
  )

  # Defining the class
  class(out) <- "ensemble_model"

  return(out)
}


#' Predictions from ensemble model
#'
#' Makes predictions using an ensemble model created from
#' \code{\link{VariableSelection}} outputs. For each observation in
#' \code{xdata}, the predictions are calculated as the average predicted values
#' obtained for that observation over the \code{K} models fitted in calibrated
#' stability selection.
#'
#' @inheritParams Ensemble
#' @param ensemble output of \code{\link{Ensemble}}.
#' @param ... additional parameters passed to \code{\link[stats]{predict}}.
#'
#' @return A matrix of predictions computed from the observations in
#'   \code{xdata}.
#'
#' @family ensemble model functions
#'
#' @examples
#' # Linear regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 50, family = "gaussian")
#' ids <- Split(data = simul$ydata, tau = c(0.8, 0.2))
#' stab <- VariableSelection(
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ]
#' )
#' ensemble <- Ensemble(
#'   stability = stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ]
#' )
#' yhat <- EnsemblePredictions(
#'   ensemble = ensemble,
#'   xdata = simul$xdata[ids[[2]], ]
#' )
#' cor(simul$ydata[ids[[2]], ], yhat)^2 # Q-squared
#'
#' # Logistic regression
#' set.seed(1)
#' simul <- SimulateRegression(n = 1000, pk = 20, family = "binomial", ev_xz = 0.9)
#' ids <- Split(data = simul$ydata, family = "binomial", tau = c(0.8, 0.2))
#' stab <- VariableSelection(
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   family = "binomial"
#' )
#' ensemble <- Ensemble(
#'   stability = stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ]
#' )
#' yhat <- EnsemblePredictions(
#'   ensemble = ensemble,
#'   xdata = simul$xdata[ids[[2]], ]
#' )
#' PlotROC(ROC(predicted = yhat, observed = simul$ydata[ids[[2]], ]))
#' yhat <- EnsemblePredictions(
#'   ensemble = ensemble,
#'   xdata = simul$xdata[ids[[2]], ],
#'   type = "response"
#' ) # predicted probabilities
#'
#' @export
EnsemblePredictions <- function(ensemble, xdata, ...) {
  # Checking family argument
  if (!ensemble$family %in% c("gaussian", "binomial")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian' or 'binomial'.")
  }

  # Checking xdata input
  if (!all(colnames(xdata) %in% colnames(ensemble$beta))) {
    stop("Arguments 'stability' and 'xdata' are not consistent. Column names in 'xdata' do not match the ones used to create 'stability'.")
  }
  xdata <- xdata[, colnames(ensemble$beta)]

  # Making predictions
  ypreds <- matrix(NA, nrow = nrow(xdata), ncol = nrow(ensemble$beta))
  for (k in 1:nrow(ensemble$beta)) {
    ypreds[, k] <- stats::predict(ensemble$models[[k]], newdata = as.data.frame(xdata), ...)
  }
  yhat <- cbind(apply(ypreds, 1, mean))
  rownames(yhat) <- rownames(xdata)

  return(yhat)
}

#' Predict method for stability selection
#'
#' Computes predicted values from the output of \code{\link{VariableSelection}}.
#'
#' @param object output of \code{\link{VariableSelection}}.
#' @param xdata predictor data (training set).
#' @param ydata outcome data (training set).
#' @param newdata optional predictor data (test set).
#' @param method character string indicating if predictions should be obtained
#'   from an \code{\link{Ensemble}} model (if \code{method="ensemble"}) or a
#'   \code{\link{Refit}}ted model (if \code{method="refit"}).
#' @param ... additional arguments passed to \code{\link[stats]{predict}}.
#'
#' @return Predicted values.
#'
#' @seealso \code{\link{Refit}}, \code{\link{Ensemble}},
#'   \code{\link{EnsemblePredictions}}
#'
#' @examples
#' \donttest{
#' ## Linear regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 50, family = "gaussian")
#'
#' # Training/test split
#' ids <- Split(data = simul$ydata, tau = c(0.8, 0.2))
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ]
#' )
#'
#' # Predictions from post stability selection estimation
#' yhat <- predict(stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   newdata = simul$xdata[ids[[2]], ],
#'   method = "refit"
#' )
#' cor(simul$ydata[ids[[2]], ], yhat)^2 # Q-squared
#'
#' # Predictions from ensemble model
#' yhat <- predict(stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   newdata = simul$xdata[ids[[2]], ],
#'   method = "ensemble"
#' )
#' cor(simul$ydata[ids[[2]], ], yhat)^2 # Q-squared
#'
#'
#' ## Logistic regression
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 500, pk = 20, family = "binomial", ev_xy = 0.9)
#'
#' # Training/test split
#' ids <- Split(data = simul$ydata, family = "binomial", tau = c(0.8, 0.2))
#'
#' # Stability selection
#' stab <- VariableSelection(
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   family = "binomial"
#' )
#'
#' # Predictions from post stability selection estimation
#' yhat <- predict(stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   newdata = simul$xdata[ids[[2]], ],
#'   method = "refit", type = "response"
#' )
#' plot(ROC(predicted = yhat, observed = simul$ydata[ids[[2]], ]))
#'
#' # Predictions from ensemble model
#' yhat <- predict(stab,
#'   xdata = simul$xdata[ids[[1]], ],
#'   ydata = simul$ydata[ids[[1]], ],
#'   newdata = simul$xdata[ids[[2]], ],
#'   method = "ensemble", type = "response"
#' )
#' plot(ROC(predicted = yhat, observed = simul$ydata[ids[[2]], ]),
#'   add = TRUE,
#'   col = "blue"
#' )
#' }
#' @export
predict.variable_selection <- function(object,
                                       xdata, ydata, newdata = NULL,
                                       method = c("ensemble", "refit"),
                                       ...) {
  # Checking inputs
  if (!object$methods$family %in% c("gaussian", "binomial", "multinomial", "cox")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian', 'binomial', 'multinomial' or 'cox'.")
  } else {
    if (method[[1]] == "ensemble") {
      if (!object$methods$family %in% c("gaussian", "binomial")) {
        method <- "refit"
        message("Predictions from ensemble models is only available for the following families for regression models: 'gaussian' or 'binomial'. Predicted values are obtained from refitting.")
      }
    }
  }

  # Using the same data if not provided
  if (is.null(newdata)) {
    newdata <- xdata
  }

  # Predictions from ensemble model
  if (method[1] == "ensemble") {
    ensemble <- Ensemble(
      stability = object,
      xdata = xdata,
      ydata = ydata
    )
    yhat <- EnsemblePredictions(
      ensemble = ensemble,
      xdata = newdata,
      ...
    )
  }

  # Predictions from refitted model
  if (method[1] == "refit") {
    refitted <- Refit(xdata = xdata, ydata = ydata, stability = object)
    if (inherits(refitted, "cv.glmnet")) {
      ids_predictors <- intersect(colnames(newdata), rownames(stats::coef(refitted)))
      yhat <- stats::predict(object = refitted, newx = as.matrix(newdata[, ids_predictors]))
    } else {
      yhat <- stats::predict(object = refitted, newdata = as.data.frame(newdata), ...)
    }
    yhat <- cbind(yhat)
  }
  return(yhat)
}

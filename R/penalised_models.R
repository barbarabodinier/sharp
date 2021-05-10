#' Penalised regression
#'
#' Runs penalised regression using implementation from
#' \code{\link[glmnet]{glmnet}}. This function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param ... additional parameters passed to \code{\link[glmnet]{glmnet}}.
#'
#' @return A list with: \item{selected}{matrix of binary selection status. Rows
#'   correspond to different model parameters. Columns correspond to
#'   predictors.} \item{beta_full}{array of model coefficients. Rows correspond
#'   to different model parameters. Columns correspond to predictors. Indices
#'   along the third dimension correspond to outcome variable(s).}
#'
#' @family underlying algorithm functions
#' @seealso \code{\link{SelectionAlgo}}, \code{\link{VariableSelection}}
#'
#' @examples
#'
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(pk = 50)
#'
#' # Running the LASSO
#' mylasso <- PenalisedRegression(
#'   xdata = simul$X, ydata = simul$Y,
#'   Lambda = c(0.1, 0.2), family = "gaussian"
#' )
#' @export
PenalisedRegression <- function(xdata, ydata, Lambda = NULL, family, ...) {
  # Making sure none of the variables has a null standard deviation
  mysd <- apply(xdata, 2, stats::sd)
  if (any(mysd == 0)) {
    for (k in which(mysd == 0)) {
      xdata[, k] <- xdata[, k] + stats::rnorm(n = nrow(xdata), sd = min(mysd[mysd != 0]) / 100)
    }
  }
  xdata <- scale(xdata)

  # Running the regression
  if (family == "multinomial") {
    mymodel <- glmnet::glmnet(x = xdata, y = ydata, lambda = Lambda, family = family, type.multinomial = "grouped", ...)
  } else {
    mymodel <- glmnet::glmnet(x = xdata, y = ydata, lambda = Lambda, family = family, ...)
  }

  if (!is.infinite(mymodel$lambda[1])) {
    # Extracting and formatting the beta coefficients
    if (!family %in% c("mgaussian", "multinomial")) {
      mybeta <- stats::coef(mymodel)
      mybeta <- t(as.matrix(mybeta))
      mybeta <- mybeta[, colnames(xdata), drop = FALSE] # removing the intercept if included

      # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
      if (any(mysd == 0)) {
        mybeta[, which(mysd == 0)] <- 0
      }

      # Preparing the outputs
      selected <- ifelse(mybeta != 0, yes = 1, no = 0)
      beta_full <- mybeta
    } else {
      if (family == "mgaussian") {
        mybeta <- array(NA,
          dim = c(length(Lambda), ncol(xdata), ncol(ydata)),
          dimnames = list(paste0("s", 0:(length(Lambda) - 1)), colnames(xdata), colnames(ydata))
        )
        for (y_id in 1:ncol(ydata)) {
          tmpbeta <- stats::coef(mymodel)[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta

          # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
          if (any(mysd == 0)) {
            mybeta[, which(mysd == 0), y_id] <- 0
          }
        }
      }
      if (family == "multinomial") {
        y_levels <- sort(unique(ydata))
        mybeta <- array(NA,
          dim = c(length(Lambda), ncol(xdata), length(y_levels)),
          dimnames = list(
            paste0("s", 0:(length(Lambda) - 1)), colnames(xdata),
            paste0("Y", y_levels)
          )
        )
        for (y_id in 1:length(y_levels)) {
          tmpbeta <- stats::coef(mymodel)[[y_id]]
          tmpbeta <- t(as.matrix(tmpbeta))
          tmpbeta <- tmpbeta[, colnames(xdata), drop = FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta), colnames(tmpbeta), y_id] <- tmpbeta

          # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
          if (any(mysd == 0)) {
            mybeta[, which(mysd == 0), y_id] <- 0
          }
        }
      }

      # Preparing the outputs
      selected <- ifelse(mybeta[, , 1, drop = FALSE] != 0, yes = 1, no = 0)
      beta_full <- mybeta
    }
  } else {
    # Returning infinite beta is the model failed
    selected <- beta_full <- Inf
  }

  return(list(selected = selected, beta_full = beta_full))
}

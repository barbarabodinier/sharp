#' @export
coef.variable_selection <- function(object, ...) {
  # Checking inputs
  if (!object$methods$family %in% c("gaussian", "binomial", "multinomial", "mgaussian", "cox")) {
    stop("This function can only be applied with the following families for regression models: 'gaussian', 'binomial', 'multinomial', 'mgaussian' or 'cox'.")
  }

  # Extracting index of calibrated parameter
  argmax_id <- ArgmaxId(stability = object)[1]

  # Extracting beta coefficients
  if (object$methods$family %in% c("gaussian", "binomial", "cox")) {
    beta <- t(object$Beta[argmax_id, , ])
    rownames(beta) <- paste0("iter", seq_len(nrow(beta)))
  }
  if (object$methods$family %in% c("multinomial", "mgaussian")) {
    tmpbeta <- object$Beta[argmax_id, , , ]
    beta <- array(NA, dim = c(dim(tmpbeta)[2], dim(tmpbeta)[1], dim(tmpbeta)[3]))
    for (k in seq_len(dim(tmpbeta)[3])) {
      beta[, , k] <- t(tmpbeta[, , k])
    }
    dimnames(beta) <- list(
      paste0("iter", seq_len(nrow(beta))),
      dimnames(object$Beta)[[2]],
      dimnames(object$Beta)[[4]]
    )
  }
  # Intercept is not included but could be obtained from Ensemble() for "gaussian" or "binomial"

  return(beta)
}


#' @export
coef.structural_model <- function(object, ...) {
  # Extracting index of calibrated parameter
  argmax_id <- ArgmaxId(stability = object)[1]

  # Extracting beta coefficients
  if (object$methods$family %in% c("gaussian", "binomial", "cox")) {
    beta <- t(object$Beta[argmax_id, , ])
    rownames(beta) <- paste0("iter", seq_len(nrow(beta)))
  }
  if (object$methods$family %in% c("multinomial", "mgaussian")) {
    tmpbeta <- object$Beta[argmax_id, , , ]
    beta <- array(NA, dim = c(dim(tmpbeta)[2], dim(tmpbeta)[1], dim(tmpbeta)[3]))
    for (k in seq_len(dim(tmpbeta)[3])) {
      beta[, , k] <- t(tmpbeta[, , k])
    }
    dimnames(beta) <- list(
      paste0("iter", seq_len(nrow(beta))),
      dimnames(object$Beta)[[2]],
      dimnames(object$Beta)[[4]]
    )
  }

  return(beta)
}

#' Classification And Regression Trees
#'
#' Runs decision trees using implementation from \code{\link[rpart]{rpart}}.
#' This function is not using stability.
#'
#' @inheritParams VariableSelection
#' @param Lambda matrix of parameters controlling the number of splits in the
#'   decision tree.
#' @param ... additional parameters passed to \code{\link[rpart]{rpart}}.
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
#' @references \insertRef{CART}{sharp}
#'
#' @examples
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(pk = 50)
#'
#'   # Running the LASSO
#'   mycart <- CART(
#'     xdata = simul$xdata,
#'     ydata = simul$ydata,
#'     family = "gaussian"
#'   )
#'   head(mycart$selected)
#' }
#' @export
CART <- function(xdata, ydata, Lambda = NULL, family, ...) {
  # Checking rpart package is installed
  CheckPackageInstalled("rpart")

  # Storing extra arguments
  extra_args <- list(...)

  # Defining default Lambda
  if (is.null(Lambda)) {
    Lambda <- cbind(seq(1, min(nrow(xdata) / 2, 100)))
  } else {
    Lambda <- cbind(Lambda)
  }

  # Defining the rpart method
  rpart_method <- switch(family,
    gaussian = "anova",
    binomial = "class",
    cox = "exp",
    poisson = "poisson"
  )

  # Writing the formula
  myformula <- stats::as.formula(paste0("ydata ~ ", paste(paste0("`", colnames(xdata), "`"), collapse = " + ")))

  # Creating the data frame
  mydata <- data.frame(ydata = ydata[, 1], xdata)

  # Initialising the parameters
  if (!"cp" %in% names(extra_args)) {
    extra_args$cp <- 0
  }
  if (!"maxdepth" %in% names(extra_args)) {
    extra_args$maxdepth <- 30
  }
  if (!any(c("minsplit", "minbucket") %in% names(extra_args))) {
    extra_args$minsplit <- 1
    extra_args$minbucket <- 1
  }

  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = rpart::rpart)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c(
    "formula", "data", "method", "model", "x", "y",
    "control", "maxcompete", "maxsurrogate", "xval", "group_x"
  )]

  # Fitting the decision tree
  mytree <- do.call(rpart::rpart, args = c(
    list(
      formula = myformula,
      data = mydata,
      method = rpart_method,
      model = FALSE,
      x = FALSE,
      y = FALSE,
      maxcompete = 0,
      maxsurrogate = 0,
      xval = 0
    ),
    tmp_extra_args
  ))

  # Pruning the decision tree
  selected <- matrix(0, nrow = nrow(Lambda), ncol = ncol(xdata))
  colnames(selected) <- colnames(xdata)
  for (i in seq_len(nrow(Lambda))) {
    n_split <- Lambda[i, 1]
    id <- max(which(mytree$cptable[, 2] <= n_split))
    mypruned <- rpart::prune(mytree, cp = mytree$cptable[id, 1])
    selected[i, unique(rownames(mypruned$splits))] <- 1
  }

  return(list(selected = selected, beta_full = selected))
}

#' Multi-block grid
#'
#' Generates a matrix of parameters controlling the sparsity of the underlying
#' selection algorithm for multi-block calibration.
#'
#' @param Lambda vector or matrix of penalty parameters.
#' @param lambda_other_blocks optional vector of penalty parameters to use for
#'   other blocks in the iterative multi-block procedure.
#'
#' @family multi-block functions
#' @seealso \code{\link{GraphicalModel}}
#'
#' @return A list with: \item{Lambda}{a matrix of (block-specific) penalty
#'   parameters. In multi-block stability selection, rows correspond to sets of
#'   penalty parameters and columns correspond to different blocks.}
#'   \item{Sequential_template}{logical matrix encoding the type of procedure
#'   for data with multiple blocks in stability selection graphical modelling.
#'   For multi-block estimation, each block is calibrated separately while
#'   others blocks are weakly penalised (\code{TRUE} only for the block
#'   currently being calibrated and \code{FALSE} for other blocks). Other
#'   approaches with joint calibration of the blocks are allowed (all entries
#'   are set to \code{TRUE}).}
#'
#' @examples
#' # Multi-block grid
#' Lambda <- matrix(
#'   c(
#'     0.8, 0.6, 0.3,
#'     0.5, 0.4, 0.2,
#'     0.7, 0.5, 0.1
#'   ),
#'   ncol = 3, byrow = TRUE
#' )
#' mygrid <- BlockLambdaGrid(Lambda, lambda_other_blocks = 0.1)
#'
#' # Multi-parameter grid (not recommended)
#' Lambda <- matrix(
#'   c(
#'     0.8, 0.6, 0.3,
#'     0.5, 0.4, 0.2,
#'     0.7, 0.5, 0.1
#'   ),
#'   ncol = 3, byrow = TRUE
#' )
#' mygrid <- BlockLambdaGrid(Lambda, lambda_other_blocks = NULL)
#' @export
BlockLambdaGrid <- function(Lambda, lambda_other_blocks = NULL) {
  if ((length(lambda_other_blocks) == 1) & (!is.vector(Lambda))) {
    lambda_other_blocks <- rep(lambda_other_blocks, ncol(Lambda))
  }
  if ((is.null(lambda_other_blocks)) & (!is.vector(Lambda))) {
    Lambda_blocks <- Lambda
    Sequential_template <- matrix(TRUE, ncol = ncol(Lambda), nrow = nrow(Lambda))
  } else {
    # Create Lambda grid matrix with nblocks columns
    if (!is.null(lambda_other_blocks)) {
      nblocks <- length(lambda_other_blocks)
    } else {
      lambda_other_blocks <- 1
      nblocks <- 1
    }
    Lambda_blocks <- NULL
    if (is.vector(Lambda)) {
      Sequential_template <- matrix(FALSE, nrow = nblocks * length(Lambda), ncol = nblocks)
    } else {
      Sequential_template <- matrix(FALSE, nrow = nblocks * nrow(Lambda), ncol = nblocks)
    }
    for (block_id in 1:nblocks) {
      if (!is.vector(Lambda)) {
        tmpLambda <- Lambda[, block_id]
      } else {
        tmpLambda <- Lambda
      }
      Lambda_blocks <- cbind(Lambda_blocks, rep(lambda_other_blocks[block_id], nblocks * length(tmpLambda)))
      Lambda_blocks[(length(tmpLambda) * (block_id - 1) + 1):(length(tmpLambda) * (block_id)), block_id] <- tmpLambda
      Sequential_template[(length(tmpLambda) * (block_id - 1) + 1):(length(tmpLambda) * (block_id)), block_id] <- TRUE
    }
  }

  return(list(Lambda = Lambda_blocks, Sequential_template = Sequential_template))
}

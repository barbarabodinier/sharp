#' sharp: Stability-enHanced Approaches using Resampling Procedures
#'
#' In stability selection and consensus clustering, resampling techniques are
#' used to enhance the reliability of the results. In this package,
#' hyper-parameters are calibrated by maximising model stability, which is
#' measured under the null hypothesis that all selection (or co-membership)
#' probabilities are identical. Functions are readily implemented for the use of
#' LASSO regression, sparse PCA, sparse (group) PLS or graphical LASSO in
#' stability selection, and hierarchical clustering, partitioning around
#' medoids, K means or Gaussian mixture models in consensus clustering.
#'
#' \tabular{ll}{ Package: \tab sharp\cr Type: \tab Package\cr Version: \tab
#' 1.4.4 \cr Date: \tab 2023-10-21 \cr License: \tab GPL (>= 3)\cr Maintainer:
#' \tab Barbara Bodinier \email{barbara.bodinier@@gmail.com}}
#'
#' @references \insertRef{OurConsensusClustering}{sharp}
#'
#' \insertRef{ourstabilityselection}{sharp}
#'
#'   \insertRef{stabilityselectionMB}{sharp}
#'
#'   \insertRef{ConsensusClustering}{sharp}
#'
#' @name sharp-package
#'
#' @importFrom Rdpack reprompt
#' @importFrom mclust mclustBIC
#'
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(5, 5, 5, 5))
#'
#' ## Regression models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' CalibrationPlot(stab)
#' summary(stab)
#' SelectedVariables(stab)
#'
#'
#' ## Graphical models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, topology = "scale-free")
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#' CalibrationPlot(stab)
#' summary(stab)
#' plot(stab)
#'
#'
#' ## PCA models
#' if (requireNamespace("elasticnet", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateComponents(pk = c(5, 3, 4))
#'   plot(simul)
#'
#'   # Stability selection
#'   stab <- BiSelection(
#'     xdata = simul$data,
#'     ncomp = 3,
#'     implementation = SparsePCA
#'   )
#'   CalibrationPlot(stab)
#'   summary(stab)
#'   SelectedVariables(stab)
#' }
#'
#'
#' ## PLS models
#' if (requireNamespace("sgPLS", quietly = TRUE)) {
#'   # Data simulation
#'   set.seed(1)
#'   simul <- SimulateRegression(n = 50, pk = c(10, 20, 30), family = "gaussian")
#'
#'   # Stability selection
#'   stab <- BiSelection(
#'     xdata = simul$xdata, ydata = simul$ydata,
#'     family = "gaussian", ncomp = 3,
#'     implementation = SparsePLS
#'   )
#'   CalibrationPlot(stab)
#'   summary(stab)
#'   plot(stab)
#' }
#'
#' par(oldpar)
#' }
NULL

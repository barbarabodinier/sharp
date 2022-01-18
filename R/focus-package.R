#' focus: Feature selectiOn and Clustering Using Stability
#'
#' Implementation of stability-enhanced regression, dimensionality reduction,
#' graphical and clustering models. These methods rely on resampling approaches
#' to estimate selection (or co-membership) probabilities. Calibration of the
#' models is done via maximisation of a stability score measuring the likelihood
#' of informative (non-uniform) selection procedure.
#'
#' \tabular{ll}{ Package: \tab focus\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2021-04-30\cr License: \tab GPL (>= 3)\cr Maintainer: \tab
#' Barbara Bodinier \email{b.bodinier@@imperial.ac.uk}}
#'
#' @references \insertRef{ourstabilityselection}{focus}
#'
#'   \insertRef{stabilityselectionSS}{focus}
#'
#'   \insertRef{stabilityselectionMB}{focus}
#'
#' @docType package
#' @name focus-package
#'
#' @importFrom Rdpack reprompt
#' @importFrom mclust mclustBIC
#'
#' @examples
#' \dontrun{
#' par(mar = c(5, 5, 5, 5))
#'
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
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateComponents(pk = c(5, 3, 4))
#' plot(simul)
#'
#' # Stability selection
#' stab <- BiSelection(
#'   xdata = simul$data,
#'   ncomp = 3,
#'   implementation = SparsePCA
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#' SelectedVariables(stab)
#'
#'
#' ## PLS models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 50, pk = c(10, 20, 30), family = "gaussian")
#'
#' # Stability selection
#' stab <- BiSelection(
#'   xdata = simul$xdata, ydata = simul$ydata,
#'   family = "gaussian", ncomp = 3,
#'   implementation = SparsePLS
#' )
#' CalibrationPlot(stab)
#' summary(stab)
#' plot(stab)
#' }
NULL

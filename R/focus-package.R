#' focus: Feature selectiOn and Clustering Using Stability
#'
#' Implementation of stability-enhanced models for variable selection in
#' regression, graphical modelling and clustering. These methods are based on
#' resampling approaches to compute selection proportions. Calibration of the
#' models is done via maximisation of a stability score measuring how unlikely
#' it is that the selection procedure is uniform.
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
#' @examples
#' \dontrun{
#'
#' ## Regression models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateRegression(n = 100, pk = 50)
#'
#' # Stability selection
#' stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
#' CalibrationPlot(stab)
#' argmax <- Argmax(stab) # calibrated parameters
#' stably_selected <- SelectedVariables(stab) # stably selected variables
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
#'   LambdaX = 1:(ncol(simul$xdata) - 1),
#'   LambdaY = 1:(ncol(simul$ydata) - 1),
#'   implementation = SparsePLS
#' )
#' print(stab$summary)
#'
#' ## Graphical models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateGraphical(n = 100, pk = 20, topology = "scale-free")
#'
#' # Stability selection
#' stab <- GraphicalModel(xdata = simul$data)
#' CalibrationPlot(stab)
#' argmax <- Argmax(stab) # calibrated parameters
#' stably_selected <- Adjacency(stab)
#' plot(Graph(stably_selected))
#'
#' ## Clustering models
#' # Data simulation
#' set.seed(1)
#' simul <- SimulateClustering(n = c(10, 10, 10), pk = 100)
#' par(mar = c(5, 5, 5, 5))
#' Heatmap(
#'   mat = cor(t(simul$data)),
#'   colours = c("navy", "white", "red"),
#'   legend_range = c(-1, 1)
#' )
#'
#' # Stability selection
#' stab <- Clustering(xdata = simul$data)
#' CalibrationPlot(stab, xlab = expression(italic(k)))
#' Clusters(stab)
#' }
NULL

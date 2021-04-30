#' focus: Feature selectiOn and Clustering Using Stability
#'
#' Implementation of stability-enhanced models for variable selection in
#' regression, network estimation and clustering. These methods are based on
#' resampling approaches to compute selection proportions. Calibration of the
#' models is done via maximisation of a stability score measuring how unlikely
#' it is that the selection procedure is uniform.
#'
#' \tabular{ll}{ Package: \tab focus\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2021-04-30\cr License: \tab GPL (>= 3)\cr }
#'
#' Very simple to use. Accepts \code{x,y} data for regression models, and
#' produces the regularization path over a grid of values for the tuning
#' parameter \code{lambda}. Only 5 functions: \code{glmnet}\cr
#' \code{predict.glmnet}\cr \code{plot.glmnet}\cr \code{print.glmnet}\cr
#' \code{coef.glmnet}
#'
#' @docType package
#' @name focus-package
#' @author Barbara Bodinier <b.bodinier@@imperial.ac.uk>
#' @examples
#' ## Regression models
#' # Data simulation
#' set.seed(1)
#' simul=SimulateRegression(n=100, pk=50)
#'
#' # Stability selection
#' stab=VariableSelection(xdata=simul$X, ydata=simul$Y)
#' CalibrationPlot(stab)
#' argmax=Argmax(stab) # calibrated parameters
#' stably_selected=SelectedVariables(stab) # stably selected variables
#'
#' ## Graphical models
#' # Data simulation
#' set.seed(1)
#' simul=SimulateGraphical(n=100, pk=20, topology="scale-free")
#'
#' # Stability selection
#' stab=GraphicalModel(data=simul$data)
#' CalibrationPlot(stab)
#' argmax=Argmax(stab) # calibrated parameters
#' stably_selected=Adjacency(stab)
#' plot(Graph(stably_selected))
NULL

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structural_equation_modelling.R
\name{StructuralModel}
\alias{StructuralModel}
\title{Stability selection in Structural Equation Modelling}
\usage{
StructuralModel(
  xdata,
  adjacency,
  residual_covariance = NULL,
  Lambda = NULL,
  pi_list = seq(0.01, 0.99, by = 0.01),
  K = 100,
  tau = 0.5,
  seed = 1,
  n_cat = NULL,
  implementation = PenalisedLinearSystem,
  resampling = "subsampling",
  cpss = FALSE,
  PFER_method = "MB",
  PFER_thr = Inf,
  FDP_thr = Inf,
  Lambda_cardinal = 100,
  optimisation = c("grid_search", "nloptr"),
  n_cores = 1,
  output_data = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{xdata}{matrix with observations as rows and variables as columns.
Column names must be defined and in line with the row and column names of
\code{adjacency}.}

\item{adjacency}{binary adjacency matrix of the Directed Acyclic Graph
(transpose of the asymmetric matrix A in Reticular Action Model notation).
The row and column names of this matrix must be defined.}

\item{residual_covariance}{binary and symmetric matrix encoding the nonzero
entries in the residual covariance matrix (symmetric matrix S in Reticular
Action Model notation). By default, this is the identity matrix (no
residual covariance).}

\item{Lambda}{matrix of parameters controlling the level of sparsity in the
underlying feature selection algorithm specified in \code{implementation}.}

\item{pi_list}{vector of thresholds in selection proportions. If
\code{n_cat=NULL} or \code{n_cat=2}, these values must be \code{>0} and
\code{<1}. If \code{n_cat=3}, these values must be \code{>0.5} and
\code{<1}.}

\item{K}{number of resampling iterations.}

\item{tau}{subsample size. Only used if \code{resampling="subsampling"} and
\code{cpss=FALSE}.}

\item{seed}{value of the seed to initialise the random number generator and
ensure reproducibility of the results (see \code{\link[base]{set.seed}}).}

\item{n_cat}{computation options for the stability score. Default is
\code{NULL} to use the score based on a z test. Other possible values are 2
or 3 to use the score based on the negative log-likelihood.}

\item{implementation}{function to use for variable selection.}

\item{resampling}{resampling approach. Possible values are:
\code{"subsampling"} for sampling without replacement of a proportion
\code{tau} of the observations, or \code{"bootstrap"} for sampling with
replacement generating a resampled dataset with as many observations as in
the full sample. Alternatively, this argument can be a function to use for
resampling. This function must use arguments named \code{data} and
\code{tau} and return the IDs of observations to be included in the
resampled dataset.}

\item{cpss}{logical indicating if complementary pair stability selection
should be done. For this, the algorithm is applied on two non-overlapping
subsets of half of the observations. A feature is considered as selected if
it is selected for both subsamples. With this method, the data is split
\code{K/2} times (\code{K} models are fitted). Only used if
\code{PFER_method="MB"}.}

\item{PFER_method}{method used to compute the upper-bound of the expected
number of False Positives (or Per Family Error Rate, PFER). If
\code{PFER_method="MB"}, the method proposed by Meinshausen and Bühlmann
(2010) is used. If \code{PFER_method="SS"}, the method proposed by Shah and
Samworth (2013) under the assumption of unimodality is used.}

\item{PFER_thr}{threshold in PFER for constrained calibration by error
control. If \code{PFER_thr=Inf} and \code{FDP_thr=Inf}, unconstrained
calibration is used (the default).}

\item{FDP_thr}{threshold in the expected proportion of falsely selected
features (or False Discovery Proportion) for constrained calibration by
error control. If \code{PFER_thr=Inf} and \code{FDP_thr=Inf}, unconstrained
calibration is used (the default).}

\item{Lambda_cardinal}{number of values in the grid of parameters controlling
the level of sparsity in the underlying algorithm. Only used if
\code{Lambda=NULL}.}

\item{optimisation}{character string indicating the type of optimisation
method. With \code{optimisation="grid_search"} (the default), all values in
\code{Lambda} are visited. Alternatively, optimisation algorithms
implemented in \code{\link[nloptr]{nloptr}} can be used with
\code{optimisation="nloptr"}. By default, we use
\code{"algorithm"="NLOPT_GN_DIRECT_L"}, \code{"xtol_abs"=0.1},
\code{"ftol_abs"=0.1} and \code{"maxeval"=Lambda_cardinal}. These values
can be changed by providing the argument \code{opts} (see
\code{\link[nloptr]{nloptr}}). For stability selection using penalised
regression, \code{optimisation="grid_search"} may be faster as it allows
for warm start.}

\item{n_cores}{number of cores to use for parallel computing (see argument
\code{workers} in \code{\link[future]{multisession}}). Using
\code{n_cores>1} is only supported with \code{optimisation="grid_search"}.}

\item{output_data}{logical indicating if the input datasets \code{xdata} and
\code{ydata} should be included in the output.}

\item{verbose}{logical indicating if a loading bar and messages should be
printed.}

\item{...}{additional parameters passed to the functions provided in
\code{implementation} or \code{resampling}.}
}
\value{
An object of class \code{variable_selection}. A list with: \item{S}{a
  matrix of the best stability scores for different parameters controlling
  the level of sparsity in the underlying algorithm.} \item{Lambda}{a matrix
  of parameters controlling the level of sparsity in the underlying
  algorithm.} \item{Q}{a matrix of the average number of selected features by
  the underlying algorithm with different parameters controlling the level of
  sparsity.} \item{Q_s}{a matrix of the calibrated number of stably selected
  features with different parameters controlling the level of sparsity.}
  \item{P}{a matrix of calibrated thresholds in selection proportions for
  different parameters controlling the level of sparsity in the underlying
  algorithm.} \item{PFER}{a matrix of upper-bounds in PFER of calibrated
  stability selection models with different parameters controlling the level
  of sparsity.} \item{FDP}{a matrix of upper-bounds in FDP of calibrated
  stability selection models with different parameters controlling the level
  of sparsity.} \item{S_2d}{a matrix of stability scores obtained with
  different combinations of parameters. Columns correspond to different
  thresholds in selection proportions.} \item{PFER_2d}{a matrix of
  upper-bounds in FDP obtained with different combinations of parameters.
  Columns correspond to different thresholds in selection proportions.}
  \item{FDP_2d}{a matrix of upper-bounds in PFER obtained with different
  combinations of parameters. Columns correspond to different thresholds in
  selection proportions.} \item{selprop}{a matrix of selection proportions.
  Columns correspond to predictors from \code{xdata}.} \item{Beta}{an array
  of model coefficients. Columns correspond to predictors from \code{xdata}.
  Indices along the third dimension correspond to different resampling
  iterations. With multivariate outcomes, indices along the fourth dimension
  correspond to outcome-specific coefficients.} \item{method}{a list with
  \code{type="variable_selection"} and values used for arguments
  \code{implementation}, \code{family}, \code{resampling}, \code{cpss} and
  \code{PFER_method}.} \item{params}{a list with values used for arguments
  \code{K}, \code{pi_list}, \code{tau}, \code{n_cat}, \code{pk}, \code{n}
  (number of observations), \code{PFER_thr}, \code{FDP_thr} and \code{seed}.
  The datasets \code{xdata} and \code{ydata} are also included if
  \code{output_data=TRUE}.} For all matrices and arrays returned, the rows
  are ordered in the same way and correspond to parameter values stored in
  \code{Lambda}.
}
\description{
Performs stability selection for Structural Equation Models. The underlying
arrow selection algorithm (e.g. regularised Structural Equation Modelling) is
run with different combinations of parameters controlling the sparsity (e.g.
penalty parameter) and thresholds in selection proportions. These two
hyper-parameters are jointly calibrated by maximisation of the stability
score.
}
\details{
In stability selection, a feature selection algorithm is fitted on
  \code{K} subsamples (or bootstrap samples) of the data with different
  parameters controlling the sparsity (\code{Lambda}). For a given (set of)
  sparsity parameter(s), the proportion out of the \code{K} models in which
  each feature is selected is calculated. Features with selection proportions
  above a threshold pi are considered stably selected. The stability
  selection model is controlled by the sparsity parameter(s) for the
  underlying algorithm, and the threshold in selection proportion:

  \eqn{V_{\lambda, \pi} = \{ j: p_{\lambda}(j) \ge \pi \} }

  In Structural Equation Modelling, "feature" refers to an arrow in the
  corresponding Directed Acyclic Graph.

  These parameters can be calibrated by maximisation of a stability score
  (see \code{\link{ConsensusScore}} if \code{n_cat=NULL} or
  \code{\link{StabilityScore}} otherwise) calculated under the null
  hypothesis of equiprobability of selection.

  It is strongly recommended to examine the calibration plot carefully to
  check that the grids of parameters \code{Lambda} and \code{pi_list} do not
  restrict the calibration to a region that would not include the global
  maximum (see \code{\link{CalibrationPlot}}). In particular, the grid
  \code{Lambda} may need to be extended when the maximum stability is
  observed on the left or right edges of the calibration heatmap. In some
  instances, multiple peaks of stability score can be observed. Simulation
  studies suggest that the peak corresponding to the largest number of
  selected features tend to give better selection performances. This is not
  necessarily the highest peak (which is automatically retained by the
  functions in this package). The user can decide to manually choose another
  peak.

  To control the expected number of False Positives (Per Family Error Rate)
  in the results, a threshold \code{PFER_thr} can be specified. The
  optimisation problem is then constrained to sets of parameters that
  generate models with an upper-bound in PFER below \code{PFER_thr} (see
  Meinshausen and Bühlmann (2010) and Shah and Samworth (2013)).

  Possible resampling procedures include defining (i) \code{K} subsamples of
  a proportion \code{tau} of the observations, (ii) \code{K} bootstrap
  samples with the full sample size (obtained with replacement), and (iii)
  \code{K/2} splits of the data in half for complementary pair stability
  selection (see arguments \code{resampling} and \code{cpss}). In
  complementary pair stability selection, a feature is considered selected at
  a given resampling iteration if it is selected in the two complementary
  subsamples.

  To ensure reproducibility of the results, the starting number of the random
  number generator is set to \code{seed}.

  For parallelisation, stability selection with different sets of parameters
  can be run on \code{n_cores} cores. Using \code{n_cores > 1} creates a
  \code{\link[future]{multisession}}. Alternatively,
  the function can be run manually with different \code{seed}s and all other
  parameters equal. The results can then be combined using
  \code{\link{Combine}}.
}
\examples{
oldpar <- par(no.readonly = TRUE)
par(mar = rep(7, 4))

\donttest{
# Data simulation
set.seed(1)
pk <- c(3, 2, 3)
simul <- SimulateStructural(
  n = 500,
  pk = pk,
  nu_between = 0.5,
  v_between = 1,
  v_sign = 1
)

# Stability selection (using glmnet)
dag <- LayeredDAG(layers = pk)
stab <- StructuralModel(
  xdata = simul$data,
  adjacency = dag
)
CalibrationPlot(stab)
LinearSystemMatrix(vect = Stable(stab), adjacency = dag)

# Stability selection (using OpenMx)
if (requireNamespace("OpenMx", quietly = TRUE)) {
  stab <- StructuralModel(
    xdata = simul$data,
    implementation = PenalisedOpenMx,
    Lambda = seq(50, 500, by = 50),
    adjacency = dag
  )
  CalibrationPlot(stab)
  OpenMxMatrix(SelectedVariables(stab), adjacency = dag)
}
}
\dontrun{
# Data simulation with latent variables
set.seed(1)
pk <- c(3, 2, 3)
simul <- SimulateStructural(
  n = 500,
  pk = pk,
  nu_between = 0.5,
  v_sign = 1,
  v_between = 1,
  n_manifest = 3,
  ev_manifest = 0.95
)

# Stability selection (using OpenMx)
if (requireNamespace("OpenMx", quietly = TRUE)) {
  dag <- LayeredDAG(layers = pk, n_manifest = 3)
  penalised <- dag
  penalised[, seq_len(ncol(simul$data))] <- 0
  stab <- StructuralModel(
    xdata = simul$data,
    implementation = PenalisedOpenMx,
    adjacency = dag,
    penalised = penalised,
    Lambda = seq(10, 100, by = 20),
    K = 10 # to increase for real use
  )
  CalibrationPlot(stab)
  ids_latent <- grep("f", colnames(dag))
  OpenMxMatrix(SelectedVariables(stab),
    adjacency = dag
  )[ids_latent, ids_latent]
}
}

par(oldpar)
}
\references{
\insertRef{JStatSoft}{sharp}

  \insertRef{ourstabilityselection}{sharp}

  \insertRef{stabilityselectionMB}{sharp}

  \insertRef{stabilityselectionSS}{sharp}

  \insertRef{RegSEM}{sharp}
}
\seealso{
\code{\link{SelectionAlgo}},
  \code{\link{Resample}}, \code{\link{StabilityScore}}

Other stability functions: 
\code{\link{BiSelection}()},
\code{\link{Clustering}()},
\code{\link{GraphicalModel}()},
\code{\link{VariableSelection}()}
}
\concept{stability functions}

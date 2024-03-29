% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda_grid.R
\name{LambdaGridRegression}
\alias{LambdaGridRegression}
\title{Grid of penalty parameters (regression model)}
\usage{
LambdaGridRegression(
  xdata,
  ydata,
  tau = 0.5,
  seed = 1,
  family = "gaussian",
  resampling = "subsampling",
  Lambda_cardinal = 100,
  check_input = TRUE,
  ...
)
}
\arguments{
\item{xdata}{matrix of predictors with observations as rows and variables as
columns.}

\item{ydata}{optional vector or matrix of outcome(s). If \code{family} is set
to \code{"binomial"} or \code{"multinomial"}, \code{ydata} can be a vector
with character/numeric values or a factor.}

\item{tau}{subsample size. Only used if \code{resampling="subsampling"} and
\code{cpss=FALSE}.}

\item{seed}{value of the seed to initialise the random number generator and
ensure reproducibility of the results (see \code{\link[base]{set.seed}}).}

\item{family}{type of regression model. This argument is defined as in
\code{\link[glmnet]{glmnet}}. Possible values include \code{"gaussian"}
(linear regression), \code{"binomial"} (logistic regression),
\code{"multinomial"} (multinomial regression), and \code{"cox"} (survival
analysis).}

\item{resampling}{resampling approach. Possible values are:
\code{"subsampling"} for sampling without replacement of a proportion
\code{tau} of the observations, or \code{"bootstrap"} for sampling with
replacement generating a resampled dataset with as many observations as in
the full sample. Alternatively, this argument can be a function to use for
resampling. This function must use arguments named \code{data} and
\code{tau} and return the IDs of observations to be included in the
resampled dataset.}

\item{Lambda_cardinal}{number of values in the grid of parameters controlling
the level of sparsity in the underlying algorithm.}

\item{check_input}{logical indicating if input values should be checked
(recommended).}

\item{...}{additional parameters passed to the functions provided in
\code{implementation} or \code{resampling}.}
}
\value{
A matrix of lambda values with one column and as many rows as
  indicated in \code{Lambda_cardinal}.
}
\description{
Generates a relevant grid of penalty parameter values for penalised
regression using the implementation in \code{\link[glmnet]{glmnet}}.
}
\examples{
# Data simulation
set.seed(1)
simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian") # simulated data

# Lambda grid for linear regression
Lambda <- LambdaGridRegression(
  xdata = simul$xdata, ydata = simul$ydata,
  family = "gaussian", Lambda_cardinal = 20
)

# Grid can be used in VariableSelection()
stab <- VariableSelection(
  xdata = simul$xdata, ydata = simul$ydata,
  family = "gaussian", Lambda = Lambda
)
print(SelectedVariables(stab))
}
\seealso{
Other lambda grid functions: 
\code{\link{LambdaGridGraphical}()},
\code{\link{LambdaSequence}()}
}
\concept{lambda grid functions}

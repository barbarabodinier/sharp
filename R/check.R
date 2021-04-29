#' Checking input parameters (regression model)
#'
#' Checks if input parameters are valid. For invalid parameters, this function
#' (i) stops the run and generates an error message, or (ii) sets the invalid
#' parameter to its default value and reports it in a warning message.
#'
#' @param xdata matrix of predictors with observations as rows and variables as
#'   columns.
#' @param ydata vector or matrix of outcome(s).
#' @param Lambda matrix of parameters controlling the underlying feature
#'   selection algorithm specified in "implementation". With
#'   implementation="glmnet", these are penalty parameters controlling the
#'   regularised model. If Lambda=NULL, \code{\link{LambdaGridRegression}} is
#'   used to define a relevant grid.
#' @param pi_list grid of values for the threshold in selection proportion. With
#'   n_cat=3, these values must be between 0.5 and 1. With n_cat=2, these values
#'   must be between 0 and 1.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param seed value of the seed to use to ensure reproducibility.
#' @param n_cat number of categories used to compute the stability score.
#'   Possible values are 2 or 3.
#' @param family type of regression model. This argument is defined as in the
#'   \code{\link{glmnet}} function from the glmnet package. Possible values
#'   include "gaussian" (linear regression), "binomial" (logistic regression),
#'   "multinomial" (multinomial regression), and "cox" (survival analysis). This
#'   argument is only used with implementation="glmnet", or with functions using
#'   the family argument in the same way (see example below).
#' @param implementation name of the function to use for variable selection.
#'   With implementation="glmnet", the function \code{\link{glmnet}} is called.
#'   Alternatively, this argument can be a character string indicating the name
#'   of a function. The function provided must use arguments called "x", "y",
#'   "lambda" and "family" and return matrices of model coefficients (see
#'   \code{\link{SelectionAlgo}}).
#' @param resampling resampling approach. Possible values are: "subsampling" for
#'   sampling without replacement of a proportion tau of the observations, or
#'   "bootstrap" for sampling with replacement generating a resampled dataset
#'   with as many observations as in the full sample. Alternatively, this
#'   argument can be a character string indicating the name of a function to use
#'   for resampling. This function must use arguments called "data" and "tau"
#'   and return IDs of observations to be included in the resampled dataset (see
#'   example in \code{\link{Resample}}).
#' @param PFER_method method used to compute the expected number of False
#'   Positives, (or Per Family Error Rate, PFER). With PFER_method="MB", the
#'   method proposed by Meinshausen and Bühlmann (2010) is used. With
#'   PFER_method="SS", the method proposed by Shah and Samworth (2013) under the
#'   assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error
#'   control. With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is
#'   used.
#' @param FDP_thr threshold in the expected proportion of falsely selected
#'   features (or False Discovery Proportion, FDP) for constrained calibration
#'   by error control. With PFER_thr=Inf and FDP_thr=Inf, unconstrained
#'   calibration is used.
#' @param Lambda_cardinal number of values in the grid.
#' @param verbose logical indicating if a message with minimum and maximum
#'   numbers of selected variables on one instance of resampled data should be
#'   printed.
#'
#' @keywords internal
CheckInputRegression=function(xdata, ydata, Lambda=NULL, pi_list=seq(0.6,0.9,by=0.01),
                              K=100, tau=0.5, seed=1, n_cat=3,
                              family="gaussian", implementation="glmnet",
                              resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                              Lambda_cardinal=100,
                              verbose=TRUE){
  # List of arguments
  myargs=c("xdata", "ydata", "Lambda", "pi_list", "K", "tau", "seed", "n_cat",
           "family", "implementation", "resampling",
           "PFER_method", "PFER_thr", "FDP_thr",
           "Lambda_cardinal", "verbose")

  # Checking the inputs (xdata and ydata)
  if (!is.null(ydata)){
    xdata=as.matrix(xdata)
    if (sum(is.na(xdata))>0){
      stop("Invalid input for argument 'xdata'. Missing values are not allowed in 'xdata'.")
    }
    if (sum(is.na(ydata))>0){
      stop("Invalid input for argument 'ydata'. Missing values are not allowed in 'ydata'.")
    }
    if ((nrow(xdata)<10)|(ncol(xdata)<=1)){
      stop("Invalid input for argument 'xdata'. Not enough data.")
    }
  }

  # Preparing xdata
  if (is.null(colnames(xdata))){
    colnames(xdata)=paste0("var", 1:ncol(xdata))
  }

  # Preparing ydata
  if (!is.null(ydata)){
    if (is.vector(ydata)|is.factor(ydata)){
      ydata=matrix(ydata, ncol=1)
    }
  }

  # Checking the inputs (xdata and ydata)
  if (!is.null(ydata)){
    if (nrow(xdata)!=nrow(ydata)){
      stop("Arguments 'xdata' and 'ydata' are not compatible. They have different numbers of observations.")
    }
  }

  # Creating dummy ydata (for resampling)
  if (is.null(ydata)){
    ydata=cbind(rep(0,nrow(xdata)))
  }

  # Naming rows of xdata and ydata
  if (is.null(rownames(xdata))&is.null(rownames(ydata))){
    rownames(xdata)=paste0("obs", 1:nrow(xdata))
    rownames(ydata)=rownames(xdata)
  } else {
    if ((is.null(rownames(xdata)))&(!is.null(rownames(ydata)))){
      rownames(xdata)=rownames(ydata)
    }
    if ((!is.null(rownames(xdata)))&(is.null(rownames(ydata)))){
      rownames(ydata)=rownames(xdata)
    }
  }

  # Re-ordering the datasets to ensure that subsamples will be the same regardless of the order of observations in the input
  ids=sort.list(rownames(xdata))
  xdata=xdata[ids,,drop=FALSE]
  ydata=ydata[ids,,drop=FALSE]

  # Further checking/preparing ydata
  if ((implementation=="glmnet")&(family=="binomial")){
    ydata=as.factor(ydata)
    if (verbose){
      print(paste0("Reference category: ", levels(ydata)[1]))
      print(paste0("Other category: ", levels(ydata)[2]))
    }
    ydata=as.numeric(ydata)-1
    ydata=matrix(ydata, ncol=1)
  }
  if ((implementation=="glmnet")&(family=="cox")){
    if ((ncol(ydata)!=2)|(length(unique(ydata[,2]))!=2)){
      stop("Invalid input for argument 'ydata'. For Cox regression using glmnet, the argument 'ydata' needs to be a matrix or data frame with two columns: the time to event and binary status.")
    }
    colnames(ydata)=c("time", "status")
    tmp=as.factor(ydata[,2])
    if (verbose){
      print(paste0("Reference category: ", levels(tmp)[1]))
      print(paste0("Other category: ", levels(tmp)[2]))
    }
    ydata[,2]=as.numeric(tmp)-1
  }
  if ((implementation=="glmnet")&(family=="multinomial")){
    if (ncol(ydata)>1){
      ydata_original=ydata
      ydata=matrix(0, nrow=nrow(ydata_original), ncol=ncol(ydata_original))
      for (j in 1:ncol(ydata)){
        tmp=as.factor(ydata_original[,j])
        if (verbose){
          print(paste0("Reference category for column ", j, ": ", levels(tmp)[1]))
          print(paste0("Other category for column ", j, ": ", levels(tmp)[2]))
        }
        ydata[,j]=(as.numeric(tmp)-1)*j
      }
      ydata=apply(ydata,1,sum)
    } else {
      ydata=as.factor(ydata)
      if (verbose){
        print(paste0("Reference category: ", levels(ydata)[1]))
        print(paste0("Other categories: ", paste(levels(ydata)[-1], collapse=", ")))
      }
      ydata=as.numeric(ydata)-1
    }
    ydata=matrix(ydata, ncol=1)
  }

  # Checking the inputs (Lambda)
  if (!is.null(Lambda)){
    if (is.matrix(Lambda)){
      Lambda_copy=Lambda
      Lambda=NULL
      for (k in 1:ncol(Lambda_copy)){
        Lambda=cbind(Lambda, as.numeric(Lambda_copy[,k]))
      }
    } else {
      Lambda=as.numeric(Lambda)
      Lambda=cbind(Lambda)
    }
    if (any(is.na(Lambda))){
      if (all(is.na(Lambda))){
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda=as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
    rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))
  }

  # Checking the inputs (pi_list)
  pi_list=sort(pi_list)
  if (n_cat==3){
    if (any(pi_list>0.5)&any(pi_list<1)){
      if ((min(pi_list)<0.5)|(max(pi_list)>1)){
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list=pi_list[which((pi_list>0.5)&(pi_list<1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list>0)&any(pi_list<1)){
      if ((min(pi_list)<0)|(max(pi_list)>1)){
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list=pi_list[which((pi_list>0)&(pi_list<1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }

  # Checking the inputs (K)
  K=as.numeric(K)
  if ((length(K)!=1)|is.na(K)){
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K=100
  }

  # Checking the inputs (tau)
  tau=as.numeric(tau)
  if ((length(tau)!=1)|is.na(tau)|(tau>=1)|(tau<=0)){
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau=0.5
  }

  # Checking the inputs (seed)
  seed=as.numeric(seed)
  if ((length(seed)!=1)|is.na(seed)){
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed=1
  }

  # Checking the inputs (n_cat)
  n_cat=as.numeric(n_cat)
  if ((length(n_cat)!=1)|is.na(n_cat)){
    warning("Invalid input for argument 'n_cat'. The argument 'seed' must be set to 2 or 3. The default value (3) was used.")
    n_cat=3
  }

  # Checking the inputs (family)
  family=as.character(family)
  if ((length(family)!=1)|is.na(family)){
    stop("Invalid input for argument 'family'. The argument 'family' must be a character string.")
  }
  if ((implementation=="glmnet")&(!family%in%c("gaussian", "binomial", "poisson", "multinomial","cox", "mgaussian"))){
    stop("Invalid input for argument 'family'. Possible values for use with the 'glmnet' implementation are: 'gaussian', 'binomial', 'poisson', 'multinomial', 'cox' or 'mgaussian'.")
  }

  # Checking the inputs (implementation)
  implementation=as.character(implementation)
  if ((length(implementation)!=1)|is.na(implementation)){
    stop("Invalid input for argument 'implementation'. This argument is the name of the function to use for variable selection, it must be a character string.")
  }

  # Checking the inputs (resampling)
  resampling=as.character(resampling)
  if ((length(resampling)!=1)|is.na(resampling)){
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }

  # Checking the inputs (PFER_method)
  PFER_method=as.character(PFER_method)
  if ((length(PFER_method)!=1)|(!PFER_method%in%c("MB","SS"))){
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  # Checking the inputs (PFER_method and resampling)
  if ((PFER_method=="SS")&(resampling=="bootstrap")){
    warning("Arguments 'resampling' and 'PFER_method' are not compatible. With 'PFER_method' set to 'SS', the resampling is done with complementary pairs of subsamples (not bootstrap).")
    resampling="subsampling"
  }

  # Checking the inputs (PFER_thr)
  PFER_thr=as.numeric(PFER_thr)
  if ((length(PFER_thr)!=1)|is.na(PFER_thr)|(PFER_thr<=0)){
    warning("Invalid input for argument 'PFER_thr'. The threshold in the upper-bound of the expected number of False Positives 'PFER_thr' must be a single positive number (or Inf). The default value (Inf) was used.")
    PFER_thr=Inf
  }

  # Checking the inputs (FDP_thr)
  FDP_thr=as.numeric(FDP_thr)
  if ((length(FDP_thr)!=1)|is.na(FDP_thr)|((!is.infinite(FDP_thr))&(FDP_thr<=0))|((!is.infinite(FDP_thr))&(FDP_thr>1))){
    warning("Invalid input for argument 'FDP_thr'. The threshold in the upper-bound of the False Discovery Proportion 'FDP_thr' must be a single number between 0 and 1 (or Inf to deactivate). The default value (Inf) was used.")
    FDP_thr=Inf
  }

  # Checking the inputs (Lambda_cardinal)
  Lambda_cardinal=as.numeric(Lambda_cardinal)
  if ((length(Lambda_cardinal)!=1)|is.na(Lambda_cardinal)|(Lambda_cardinal<1)){
    warning("Invalid input for argument 'Lambda_cardinal'. The argument 'Lambda_cardinal' must be a single positive number. The default value (100) was used.")
    Lambda_cardinal=100
  }

  # Checking the inputs (verbose)
  verbose=as.logical(verbose)
  if ((length(verbose)!=1)|is.na(verbose)){
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose=TRUE
  }

  # Assigning checked values to the parent function
  for (i in 1:length(myargs)){
    if (!is.null(get(myargs[i]))){
      assign(myargs[i], get(myargs[i]), envir=parent.frame(n=1))
    }
  }
}


#' Checking input parameters (graphical model)
#'
#' Checks if input parameters are valid.
#' For invalid parameters, this function (i) stops the run and generates
#' an error message, or (ii) sets the invalid parameter to its default value
#' and reports it in a warning message.
#'
#' @param data matrix with observations as rows and variables as columns.
#' @param pk vector encoding the grouping structure.
#' Only used for multi-block stability selection.
#' For this, the variables in data have to be ordered
#' by group and argument "pk" has to be a vector
#' indicating the number of variables
#' in each of the groups (see example below).
#' If pk=NULL, single-block stability selection is performed.
#' @param Lambda matrix of parameters controlling the underlying
#' feature selection algorithm specified in "implementation".
#' With implementation="glassoFast", these are penalty parameters
#' controlling the regularised model.
#' If Lambda=NULL, \code{\link{LambdaGridGraphical}} is used to define
#' a relevant grid.
#' For multi-block calibration (i.e. when argument "pk" is a vector),
#' Lambda can be a vector, to use the procedure from Equation (5) (recommended),
#' or a matrix with as many columns as there are entries in "pk",
#' to use the procedure from Equation (4) (see details and examples below).
#' @param lambda_other_blocks optional vector of (penalty) parameters
#' to use for other blocks
#' in the iterative multi-block procedure (see example below).
#' Only used for multi-block graphical models, i.e. when pk is a vector.
#' @param pi_list grid of values for the threshold in selection proportion.
#' With n_cat=3, these values must be between 0.5 and 1.
#' With n_cat=2, these values must be between 0 and 1.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param seed value of the seed to use to ensure reproducibility.
#' @param n_cat number of categories used to compute the stability score.
#' Possible values are 2 or 3.
#' @param implementation name of the function to use for graphical modelling.
#' With implementation="glassoFast", the function \code{\link{glassoFast}}
#' is used for regularised estimation of a conditional independence graph.
#' Alternatively, this argument can be a character string indicating the name of a function.
#' The function provided must use arguments called "x", "lambda" and "scale"
#' and return a binary and symmetric adjacency matrix (see \code{\link{GraphicalAlgo}}).
#' @param start character string indicating if the algorithm should be
#' initialised at the estimated (inverse) covariance with previous
#' penalty parameters (start="warm") or not (start="cold").
#' Using start="warm" can speed-up the computations.
#' Only used for implementation="glassoFast" (see argument "start"
#' in \code{\link{glassoFast}}).
#' @param scale logical indicating if the correlation (if scale=TRUE)
#' or covariance (if scale=FALSE) matrix should be used as input
#' for the graphical LASSO. If implementation is not set to "glassoFast",
#' this argument must be used as input of the function provided instead.
#' @param resampling resampling approach. Possible values are: "subsampling"
#' for sampling without replacement of a proportion tau of the observations, or
#' "bootstrap" for sampling with replacement generating a resampled dataset with
#' as many observations as in the full sample. Alternatively, this argument can be
#' a character string indicating the name of a function to use for resampling.
#' This function must use arguments called "data" and "tau" and return
#' IDs of observations to be included in the resampled dataset
#' (see example in \code{\link{Resample}}).
#' @param PFER_method method used to compute the expected number of False Positives,
#' (or Per Family Error Rate, PFER). With PFER_method="MB", the method
#' proposed by Meinshausen and Bühlmann (2010) is used. With PFER_method="SS",
#' the method proposed by Shah and Samworth (2013) under the assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error control.
#' With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is used.
#' The grid is defined such that the estimated graph does not generate an upper-bound in
#' PFER above PFER_thr.
#' @param FDP_thr threshold in the expected proportion of falsely selected edges
#' (or False Discovery Proportion, FDP)
#' for constrained calibration by error control.
#' With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is used.
#' If FDP_thr is not infinite, the grid is defined such that the estimated graph
#' does not generate an upper-bound in PFER above the number of node pairs.
#' @param Lambda_cardinal number of values in the grid.
#' @param lambda_max maximum value in the grid. With lambda_max=NULL,
#' the maximum value is set to the maximum covariance in absolute value.
#' @param lambda_path_factor multiplicative factor used to define the minimum value in the grid.
#' @param max_density threshold on the density. The grid is defined such that the density
#' of the estimated graph does not exceed max_density.
#' @param verbose logical indicating if a message with minimum and maximum
#' numbers of selected variables on one instance of resampled data should be printed.
#'
#' @keywords internal
CheckInputGraphical=function(data, pk=NULL, Lambda=NULL, lambda_other_blocks=NULL,
                           pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, n_cat=3,
                           implementation="glassoFast", start="cold", scale=TRUE,
                           resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                           Lambda_cardinal=50, lambda_max=NULL, lambda_path_factor=0.0001, max_density=0.3,
                           verbose=TRUE){
  # List of arguments
  myargs=c("data", "pk", "Lambda", "lambda_other_blocks",
           "pi_list", "K", "tau", "seed", "n_cat",
           "implementation", "start", "scale",
           "resampling", "PFER_method", "PFER_thr", "FDP_thr",
           "Lambda_cardinal",
           "lambda_path_factor", "max_density",
           "verbose")

  # Checking the inputs (data)
  data=as.matrix(data)
  if (sum(is.na(data))>0){
    stop("Invalid input for argument 'data'. Missing values are not allowed in 'data'.")
  }
  if ((nrow(data)<10)|(ncol(data)<=1)){
    stop("Invalid input for argument 'data'. Not enough data.")
  }

  # Checking the inputs (pk)
  if (!is.null(pk)){
    pk=as.numeric(pk)
    if (sum(pk)!=ncol(data)){
      stop("Invalid input for argument 'pk'. The number of variables per group 'pk' must sum to the number of columns in 'data'.")
    }
  } else {
    pk=ncol(data)
  }

  # Checking the inputs (pi_list)
  pi_list=sort(pi_list)
  if (n_cat==3){
    if (any(pi_list>0.5)&any(pi_list<1)){
      if ((min(pi_list)<0.5)|(max(pi_list)>1)){
        warning("The values in 'pi_list' must be between 0.5 and 1. All other values were discarded.")
        pi_list=pi_list[which((pi_list>0.5)&(pi_list<1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0.5 and lower than 1. To consider thresholds below 0.5, argument 'n_cat' must be set to 2.")
    }
  } else {
    if (any(pi_list>0)&any(pi_list<1)){
      if ((min(pi_list)<0)|(max(pi_list)>1)){
        warning("The values in 'pi_list' must be between 0 and 1. All other values were discarded.")
        pi_list=pi_list[which((pi_list>0)&(pi_list<1))]
      }
    } else {
      stop("Invalid input for argument 'pi_list'. The values in the vector must be greater than 0 and lower than 1.")
    }
  }

  # Checking the inputs (K)
  K=as.numeric(K)
  if ((length(K)!=1)|is.na(K)){
    warning("Invalid input for argument 'K'. The number of resampling iterations 'K' must be a single number.")
    K=100
  }

  # Checking the inputs (tau)
  tau=as.numeric(tau)
  if ((length(tau)!=1)|is.na(tau)|(tau>=1)|(tau<=0)){
    warning("Invalid input for argument 'tau'. The subsample size 'tau' must be a number between 0 and 1. The default value (0.5) was used.")
    tau=0.5
  }

  # Checking the inputs (seed)
  seed=as.numeric(seed)
  if ((length(seed)!=1)|is.na(seed)){
    warning("Invalid input for argument 'seed'. The argument 'seed' must be a single number. The default value (1) was used.")
    seed=1
  }

  # Checking the inputs (n_cat)
  n_cat=as.numeric(n_cat)
  if ((length(n_cat)!=1)|is.na(n_cat)){
    warning("Invalid input for argument 'n_cat'. The argument 'seed' must be set to 2 or 3. The default value (3) was used.")
    n_cat=3
  }

  # Checking the inputs (implementation)
  implementation=as.character(implementation)
  if ((length(implementation)!=1)|is.na(implementation)){
    stop("Invalid input for argument 'implementation'. This argument is the name of the function to use for network estimation, it must be a character string.")
  }

  # Checking the inputs (start)
  start=as.character(start)
  if ((length(start)!=1)|is.na(start)|(!start%in%c("cold","warm"))){
    warning("Invalid input for argument 'start'. The argument must be 'cold' or 'warm'. The default value (cold) was used.")
  }

  # Checking the inputs (scale)
  scale=as.logical(scale)
  if ((length(scale)!=1)|is.na(scale)){
    stop("Invalid input for argument 'scale'. The argument 'scale' must be logical (TRUE or FALSE).")
  }

  # Checking the inputs (resampling)
  resampling=as.character(resampling)
  if ((length(resampling)!=1)|is.na(resampling)){
    stop("Invalid input for argument 'resampling'. The argument 'resampling' must be a character string. Possible values are: 'subsampling', 'bootstrap' or the name of a function.")
  }

  # Checking the inputs (PFER_method)
  PFER_method=as.character(PFER_method)
  if ((length(PFER_method)!=1)|(!PFER_method%in%c("MB","SS"))){
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  # Checking the inputs (PFER_method and resampling)
  if ((PFER_method=="SS")&(resampling=="bootstrap")){
    warning("Arguments 'resampling' and 'PFER_method' are not compatible. With 'PFER_method' set to 'SS', the resampling is done with complementary pairs of subsamples (not bootstrap).")
    resampling="subsampling"
  }

  # Checking the inputs (lambda_max)
  if (!is.null(lambda_max)){
    lambda_max=as.numeric(lambda_max)
    if ((length(lambda_max)!=1)|is.na(lambda_max)|(lambda_max<=0)){
      warning("Invalid input for argument 'lambda_max'. The argument 'lambda_max' must be a single positive number. The default value (NULL) was used.")
      lambda_max=NULL
    }
  }

  # Checking the inputs (lambda_path_factor)
  lambda_path_factor=as.numeric(lambda_path_factor)
  if ((length(lambda_path_factor)!=1)|is.na(lambda_path_factor)|(lambda_path_factor<=0)|(lambda_path_factor>=1)){
    warning("Invalid input for argument 'lambda_path_factor'. The argument 'lambda_path_factor' must be a single number between 0 and 1. The default value (0.0001) was used.")
    lambda_path_factor=0.0001
  }

  # Checking the inputs (max_density)
  max_density=as.numeric(max_density)
  if ((length(max_density)!=1)|is.na(max_density)|(max_density<=0)|(max_density>1)){
    warning("Invalid input for argument 'max_density'. The argument 'max_density' must be a single number between 0 and 1. The default value (0.3) was used.")
    max_density=0.3
  }

  # Checking the inputs (Lambda_cardinal)
  Lambda_cardinal=as.numeric(Lambda_cardinal)
  if ((length(Lambda_cardinal)!=1)|is.na(Lambda_cardinal)|(Lambda_cardinal<1)){
    warning("Invalid input for argument 'Lambda_cardinal'. The argument 'Lambda_cardinal' must be a single positive number. The default value (50) was used.")
    Lambda_cardinal=50
  }

  # Create matrix with block indices
  bigblocks=BlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)

  # Checking the inputs (lambda_other_blocks in single-block analyses)
  if (!is.null(lambda_other_blocks)){
    if ((length(pk)==1)){
      warning("Unused argument 'lambda_other_blocks'. This argument is specific to multi-block analyses (i.e. with multiple entries in argument 'pk').")
      lambda_other_blocks=NULL
    }
    if (length(lambda_other_blocks)==1){
      lambda_other_blocks=rep(lambda_other_blocks, nblocks)
    } else {
      if (length(lambda_other_blocks)!=nblocks){
        stop(paste0("Invalid input for argument 'lambda_other_blocks'. This argument must be a vector with as many entries as there are blocks in the data (i.e. ",
                    nblocks, " entries in this case)."))
      }
    }
  }

  # Checking the inputs (verbose)
  verbose=as.logical(verbose)
  if ((length(verbose)!=1)|is.na(verbose)){
    warning("Invalid input for argument 'verbose'. The argument 'verbose' must be logical (TRUE or FALSE). The default value (TRUE) was used.")
    verbose=TRUE
  }

  # Checking the inputs (Lambda)
  if (!is.null(Lambda)){
    if (is.matrix(Lambda)){
      if ((ncol(Lambda)!=nblocks)&(ncol(Lambda)!=1)){
        stop(paste0("Invalid input for argument 'Lambda'. The argument 'Lambda' must be a matrix as many columns as blocks (N=",nblocks,")."))
      }
      if (ncol(Lambda)==1){
        Lambda=as.numeric(as.vector(Lambda))
      } else {
        Lambda_copy=Lambda
        Lambda=NULL
        for (k in 1:ncol(Lambda_copy)){
          Lambda=cbind(Lambda, as.numeric(Lambda_copy[,k]))
        }
      }
    } else {
      Lambda=as.numeric(Lambda)
    }
    if (any(is.na(Lambda))){
      if (all(is.na(Lambda))){
        stop("Invalid input for argument 'Lambda'. The input only contains missing values.")
      } else {
        Lambda=as.matrix(stats::na.exclude(Lambda))
        warning("Invalid input for argument 'Lambda'. The input contains missing values. These have been excluded.")
      }
    }
  }

  # Checking the inputs (PFER_thr)
  PFER_thr=as.numeric(PFER_thr)
  if ((!length(PFER_thr)%in%c(1,nblocks))|is.na(PFER_thr)|(PFER_thr<=0)){
    warning("Invalid input for argument 'PFER_thr'. The threshold in the upper-bound of the expected number of False Positives 'PFER_thr' must be a vector with positive numbers (or Inf). The default value (Inf) was used.")
    PFER_thr=Inf
  }

  # Checking the inputs (FDP_thr)
  FDP_thr=as.numeric(FDP_thr)
  if (length(pk)==1){
    if ((!length(PFER_thr)%in%c(1,nblocks))|is.na(FDP_thr)|((!is.infinite(FDP_thr))&(FDP_thr<=0))|((!is.infinite(FDP_thr))&(FDP_thr>1))){
      warning("Invalid input for argument 'FDP_thr'. The threshold in the upper-bound of the False Discovery Proportion 'FDP_thr' must be a vector with numbers between 0 and 1 (or Inf to deactivate). The default value (Inf) was used.")
      FDP_thr=Inf
    }
  }

  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    }
  }

  # Assigning checked values to the parent function
  for (i in 1:length(myargs)){
    assign(myargs[i], get(myargs[i]), envir=parent.frame(n=1))
  }

  # Assigning extra objects to the parent function
  myextra=c("bigblocks", "bigblocks_vect", "blocks", "N_blocks", "nblocks", "PFER_thr_blocks", "FDP_thr_blocks")
  for (i in 1:length(myextra)){
    assign(myextra[i], get(myextra[i]), envir=parent.frame(n=1))
  }
}



#' Grid of penalty parameters (regression model)
#'
#' Generates a relevant grid of penalty parameter values for penalised
#' regression.
#'
#' @param xdata matrix of predictors with observations as rows and variables as
#'   columns.
#' @param ydata vector or matrix of outcome(s).
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param seed used in set.seed() to ensure reproducibility of the analyses.
#' @param family type of regression model. This argument is defined as in the
#'   \code{\link[glmnet]{glmnet}} function from the glmnet package. Possible
#'   values include "gaussian" (linear regression), "binomial" (logistic
#'   regression), "multinomial" (multinomial regression), and "cox" (survival
#'   analysis). This argument is only used with implementation="glmnet", or with
#'   functions using the family argument in the same way.
#' @param implementation name of the function to use for definition of the grid
#'   of lambda values. With implementation="glmnet", the function
#'   \code{\link[glmnet]{glmnet}} is used to define the path of lambda values.
#'   Alternatively, this argument can be a character string indicating the name
#'   of a function. The function provided must use arguments called "x", "y" and
#'   "family" and return a list in which the entry named "lambda" contains a
#'   vector of lambda values.
#' @param resampling resampling approach. Possible values are: "subsampling" for
#'   sampling without replacement of a proportion tau of the observations, or
#'   "bootstrap" for sampling with replacement generating a resampled dataset
#'   with as many observations as in the full sample. Alternatively, this
#'   argument can be a character string indicating the name of a function to use
#'   for resampling. This function must use arguments called "data" and "tau"
#'   and return IDs of observations to be included in the resampled dataset (see
#'   example in \code{\link{Resample}}).
#' @param Lambda_cardinal number of values in the grid.
#' @param verbose logical indicating if a message with minimum and maximum
#'   numbers of selected variables on one instance of resampled data should be
#'   printed.
#' @param ... additional parameters passed to the functions provided in
#'   "implementation" or "resampling".
#'
#' @return a matrix of lambda values with one column and as many rows as
#'   indicated in "Lambda_cardinal".
#'
#' @family lambda grid functions
#'
#' @examples
#' # Lambda grid for linear regression
#' simul=SimulateRegression(n=100, pk=20, family="gaussian") # simulated data
#' Lambda=LambdaGridRegression(xdata=simul$X, ydata=simul$Y,
#' family="gaussian", Lambda_cardinal=20)
#'
#' # Grid can be used in VariableSelection()
#' out=VariableSelection(xdata=simul$X, ydata=simul$Y,
#' family="gaussian", Lambda=Lambda)
#' SelectedVariables(out)
#'
#' # For use with gglasso (group LASSO)
#' require(gglasso)
#' ManualGridGroupLasso=function(x, y, family, ...){
#' if (family=="gaussian"){
#' return(cv.gglasso(x=x, y=y, pred.loss="L1", ...))
#' }
#' }
#' Lambda=LambdaGridRegression(xdata=simul$X, ydata=simul$Y,
#' family="gaussian", Lambda_cardinal=20,
#' implementation="ManualGridGroupLasso",
#' group=rep(1:4, each=5))
#'
#' @export
LambdaGridRegression=function(xdata, ydata, tau=0.5, seed=1,
                              family="gaussian", implementation="glmnet",
                              resampling="subsampling",
                              Lambda_cardinal=100,
                              verbose=TRUE, ...){
  # Object preparation, error and warning messages
  Lambda=NULL
  pi_list=seq(0.6,0.9,by=0.01)
  K=100
  n_cat=3
  PFER_method="MB"
  PFER_thr=Inf
  FDP_thr=Inf
  CheckInputRegression(xdata=xdata, ydata=ydata, Lambda=Lambda, pi_list=pi_list,
                       K=K, tau=tau, seed=seed, n_cat=n_cat,
                       family=family, implementation=implementation,
                       resampling=resampling, PFER_method=PFER_method,
                       PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                       Lambda_cardinal=Lambda_cardinal,
                       verbose=verbose)
  rm(n_cat)
  rm(Lambda)
  rm(pi_list)
  rm(K)

  # Taking one subsample/boostrap sample of the data
  set.seed(1) # To keep to allow for reproducible parallelisation
  s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

  # Getting upperbound of Lambda
  set.seed(1) # To keep to allow for reproducible parallelisation
  if (implementation=="glmnet"){
    mycv=glmnet::glmnet(x=xdata[s,], y=ydata[s,], family=family, ...)
  } else {
    # Applying user-defined function for variable selection
    mycv=do.call(get(implementation), args=list(x=xdata[s,], y=ydata[s,], family=family, ...))
  }

  # Printing messages
  if (verbose){
    if (implementation=="glmnet"){
      print(paste("Minimum number of selected variables:", min(mycv$nzero)))
      print(paste("Maximum number of selected variables:", max(mycv$nzero)))
    }
  }

  # Creating a grid of lambda values from min and max
  Lambda=cbind(LambdaSequence(lmax=max(mycv$lambda), lmin=min(mycv$lambda), cardinal=Lambda_cardinal))
  Lambda=as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))

  return(Lambda)
}


#' Grid of penalty parameters (graphical model)
#'
#' Generates a relevant grid of penalty parameter values for penalised graphical
#' models.
#'
#' @param data matrix with observations as rows and variables as columns.
#' @param pk vector encoding the grouping structure. Only used for multi-block
#'   stability selection. For this, the variables in data have to be ordered by
#'   group and argument "pk" has to be a vector indicating the number of
#'   variables in each of the groups. If pk=NULL, single-block stability
#'   selection is performed.
#' @param lambda_other_blocks vector of penalty parameters to use for other
#'   blocks in the iterative multi-block procedure. Only used for multi-block
#'   graphical models, i.e. when pk is not set to NULL.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param n_cat number of categories used to compute the stability score.
#'   Possible values are 2 or 3.
#' @param implementation name of the function to use for definition of the grid
#'   of lambda values. With implementation="glassoFast", the function
#'   \code{\link[glassoFast]{glassoFast}} is called and iteratively applied on
#'   possible penalty values until the constraint are verified, i.e. that the
#'   expected density is below the value given in "max_density", that the
#'   expected PFER is below the value given in "PFER_thr" or that the expected
#'   PFER is below the number of selected edges if "FDP_thr" is not set to Inf.
#'   Alternatively, this argument can be a character string indicating the name
#'   of a function. The function provided must use arguments called "x",
#'   "lambda" and "scale" and return a binary and symmetric adjacency matrix.
#' @param start character string indicating if the algorithm should be
#'   initialised at the estimated (inverse) covariance with previous penalty
#'   parameters (start="warm") or not (start="cold"). Using start="warm" can
#'   speed-up the computations. Only used for implementation="glassoFast" (see
#'   argument "start" in \code{\link[glassoFast]{glassoFast}}).
#' @param scale logical indicating if the correlation (if scale=TRUE) or
#'   covariance (if scale=FALSE) matrix should be used as input for the
#'   graphical LASSO. If implementation is not set to "glassoFast", this
#'   argument must be used as input of the function provided instead.
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
#'   method proposed by Meinshausen and Buhlmann (2010) is used. With
#'   PFER_method="SS", the method proposed by Shah and Samworth (2013) under the
#'   assumption of unimodality is used.
#' @param PFER_thr threshold in PFER for constrained calibration by error
#'   control. With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is
#'   used. The grid is defined such that the estimated graph does not generate
#'   an upper-bound in PFER above PFER_thr.
#' @param FDP_thr threshold in the expected proportion of falsely selected edges
#'   (or False Discovery Proportion, FDP) for constrained calibration by error
#'   control. With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is
#'   used. If FDP_thr is not infinite, the grid is defined such that the
#'   estimated graph does not generate an upper-bound in PFER above the number
#'   of node pairs.
#' @param Lambda_cardinal number of values in the grid.
#' @param lambda_max maximum value in the grid. With lambda_max=NULL, the
#'   maximum value is set to the maximum covariance in absolute value.
#' @param lambda_path_factor multiplicative factor used to define the minimum
#'   value in the grid. The grid is defined iteratively by multiplying the
#'   smallest value at current iteration by \code{lambda_path_factor} until one
#'   of the \code{max_density} or \code{PFER_thr} stopping criteria is met.
#' @param max_density threshold on the density. The grid is defined such that
#'   the density of the estimated graph does not exceed \code{max_density}.
#' @param verbose logical indicating if a message with minimum and maximum
#'   numbers of selected variables on one instance of resampled data should be
#'   printed.
#' @param ... additional parameters passed to the functions provided in
#'   "implementation" or "resampling".
#'
#' @return a matrix of lambda values with as many columns as there are entries
#'   in "pk" and as many rows as indicated in "Lambda_cardinal".
#'
#' @family lambda grid functions
#'
#' @references \insertRef{stabilityselectionMB}{focus}
#'
#'   \insertRef{stabilityselectionSS}{focus}
#'
#' @examples
#' # Single-block simulation
#' set.seed(1)
#' simul=SimulateGraphical()
#'
#' # Generating grid of 10 values
#' Lambda=LambdaGridGraphical(data=simul$data, Lambda_cardinal=10)
#'
#' # Ensuring PFER < 5
#' Lambda=LambdaGridGraphical(data=simul$data, Lambda_cardinal=10, PFER_thr=5)
#'
#' # Multi-block simulation
#' set.seed(1)
#' simul=SimulateGraphical(pk=c(10,10))
#'
#' # Multi-block grid
#' Lambda=LambdaGridGraphical(data=simul$data, pk=c(10,10), Lambda_cardinal=10)
#'
#' # Denser neighbouring blocks
#' Lambda=LambdaGridGraphical(data=simul$data, pk=c(10,10),
#' Lambda_cardinal=10, lambda_other_blocks=0)
#'
#' # Using different neighbour penalties
#' Lambda=LambdaGridGraphical(data=simul$data, pk=c(10,10),
#' Lambda_cardinal=10, lambda_other_blocks=c(0.1, 0, 0.1))
#' stab=GraphicalModel(data=simul$data, pk=c(10,10),
#' Lambda=Lambda, lambda_other_blocks=c(0.1,0,0.1))
#' stab$Lambda
#'
#' # Visiting from empty to full graphs with max_density=1
#' Lambda=LambdaGridGraphical(data=simul$data, pk=c(10,10),
#' Lambda_cardinal=10, max_density=1)
#' bigblocks=BlockMatrix(pk=c(10,10))
#' bigblocks_vect=bigblocks[upper.tri(bigblocks)]
#' N_blocks=unname(table(bigblocks_vect))
#' N_blocks # max number of edges per block
#' stab=GraphicalModel(data=simul$data, pk=c(10,10),Lambda=Lambda)
#' apply(stab$Q,2,max,na.rm=TRUE) # max average number of edges from underlying algo
#'
#' @export
LambdaGridGraphical=function(data, pk=NULL, lambda_other_blocks=0.1, K=100, tau=0.5, n_cat=3,
                           implementation="glassoFast", start="cold", scale=TRUE,
                           resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                           Lambda_cardinal=50, lambda_max=NULL, lambda_path_factor=0.001, max_density=0.5,
                           verbose=TRUE, ...){
  # K to keep for PFER computations with PFER_method set to "SS"
  # Error and warning messages
  bigblocks=bigblocks_vect=blocks=N_blocks=nblocks=PFER_thr_blocks=FDP_thr_blocks=NULL
  Lambda=NULL
  seed=1 # To keep to allow for reproducible parallelisation
  pi_list=0.75 # only used for screening
  CheckInputGraphical(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=lambda_other_blocks,
                    pi_list=pi_list, K=K, tau=tau, seed=seed, n_cat=n_cat,
                    implementation=implementation, start=start, scale=scale,
                    resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                    Lambda_cardinal=Lambda_cardinal,
                    lambda_max=lambda_max, lambda_path_factor=lambda_path_factor, max_density=max_density,
                    verbose=verbose)
  rm(Lambda)

  # Preparing lambda_dense
  ldense=lambda_other_blocks
  p=sum(pk)
  N=p*(p-1)/2

  # Printing message
  if (verbose){
    if (all(!is.infinite(PFER_thr_blocks))){
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))){
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }

  # Making sure none of the variables has a null standard deviation
  mysd=apply(data,2,stats::sd)
  if (any(mysd==0)){
    for (k in which(mysd==0)){
      data[,k]=data[,k]+stats::rnorm(n=nrow(data), sd=min(mysd[mysd!=0])/100)
    }
  }

  # Get upperbound of Lambda
  if (scale){
    mycov=stats::cor(data)
  } else {
    mycov=stats::cov(data)
  }

  # Theoretical starting point for lambda
  if (is.null(lambda_max)){
    diag(mycov)=0
    lambda_max=max(abs(mycov))
  }
  lmin=lambda_max
  lmin=rep(lmin, nblocks)

  # Identifying the constraint
  if (all(is.infinite(PFER_thr))){
    type_opt_problem="unconstrained"
    if (!all(is.infinite(FDP_thr))){
      type_opt_problem="constrained_PFER"
      PFER_thr=N # very loose stopping criterion for constraint on the FDP
    }
  } else {
    type_opt_problem="constrained_PFER"
  }

  if (type_opt_problem=="unconstrained"){
    max_q=rep(0,nblocks)
    redo=TRUE
    done=rep(0, nblocks)
    while (redo){
      lmin=lmin*lambda_path_factor
      Lambda=NULL
      for (b in 1:nblocks){
        Lambda=cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal=Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[l,,drop=FALSE]
        myscreen=SerialGraphical(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1,
                                       tau=tau, seed=seed, n_cat=n_cat,
                                       implementation=implementation, start=start, scale=scale,
                                       resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                                       verbose=FALSE, ...) # Only 1 iteration to get the Q

        if (l<nrow(Lambda)){
          # Updating the smallest lambda if the density of the block is still below max_density
          for (b in 1:nblocks){
            lmin[b]=ifelse((myscreen$Q[b,b]<(max_density*N_blocks)[b])&(done[b]==0),
                           yes=Lambda[l+1,b], no=lmin[b])
            done[b]=ifelse(myscreen$Q[b,b]>=(max_density*N_blocks)[b], yes=1, no=0)
          }
        }

        # Increment if max_density is not yet reached
        Q_block_iteration=NULL
        for (b in 1:nblocks){
          Q_block_iteration=c(Q_block_iteration, myscreen$Q[b,b])
        }

        if (any(Q_block_iteration<max_density*N_blocks)){
          l=l+1
        } else {
          l=nrow(Lambda) # stopping current while loop
          redo=FALSE # stopping overarching while loop
        }
      }
    }
  }

  if (type_opt_problem=="constrained_PFER"){
    max_q=rep(0,nblocks)
    redo=TRUE
    done=rep(0, nblocks)
    while (redo){
      lmin=lmin*lambda_path_factor
      Lambda=NULL
      for (b in 1:nblocks){
        Lambda=cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal=Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[l,,drop=FALSE]
        myscreen=SerialGraphical(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1,
                                       tau=tau, seed=seed, n_cat=n_cat,
                                       resampling=resampling, scale=scale,
                                       implementation=implementation, start=start, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                                       verbose=FALSE, ...) # Only 1 iteration to get the Q

        # Compute PFER
        PFER_l=rep(NA, nblocks)
        for (b in 1:nblocks){
          mytmplist=NULL
          for (j in 1:length(pi_list)){
            pi=pi_list[j]
            mytmplist=c(mytmplist, PFER(q=myscreen$Q[b,b], pi=pi, N=N_blocks[b], K=K, PFER_method=PFER_method))
          }
          PFER_l[b]=min(mytmplist)
        }

        if (l<nrow(Lambda)){
          # Updating the smallest lambda if the PFER of the block is still below the threshold (with some margin)
          lmin=ifelse((PFER_l<=(PFER_thr_blocks*1.2+1))&(done==0), yes=Lambda[l+1,], no=lmin)
          done=ifelse(PFER_l>(PFER_thr_blocks*1.2+1), yes=1, no=0)
        }

        # Increment if PFER or max_density are not yet reached
        Q_block_iteration=NULL
        for (b in 1:nblocks){
          Q_block_iteration=c(Q_block_iteration, myscreen$Q[b,b])
        }

        if (any(PFER_l<=(PFER_thr_blocks*1.2+1))&(any(Q_block_iteration<max_density*N_blocks))){
          l=l+1
        } else {
          l=nrow(Lambda) # stopping current while loop
          redo=FALSE # stopping overarching while loop
        }
      }
    }
  }

  # Prepare final lambda path for each block
  Lambda=NULL
  for (b in 1:nblocks){
    Lambda=cbind(Lambda, LambdaSequence(lambda_max, lmin[b], cardinal=Lambda_cardinal))
  }
  Lambda=as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))

  return(Lambda)
}


#' Sequence of penalty parameters
#'
#' Generates a sequence of penalty parameters from extreme values and the
#' required number of elements. The sequence is defined on the log-scale.
#'
#' @param lmax maximum value in the grid.
#' @param lmin minimum value in the grid.
#' @param cardinal number of values in the grid.
#'
#' @return a vector with values between "lmin" and "lmax" and as many values as
#'   indicated by "cardinal".
#'
#' @family lambda grid functions
#'
#' @export
LambdaSequence=function(lmax, lmin, cardinal=100){
  return(exp(seq(log(lmax), log(lmin), length.out=cardinal)))
  # return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}


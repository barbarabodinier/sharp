#' Stability selection in regression
#'
#' Runs stability selection regression models with different combinations of
#' parameters controlling the sparsity of the underlying selection algorithm
#' (e.g. penalty parameter for regularised models) and thresholds in selection
#' proportions. These two parameters are jointly calibrated by maximising the
#' stability score of the model (possibly under a constraint on the expected
#' number of falsely stably selected features).
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
#'   \code{\link[glmnet]{glmnet}} function from the glmnet package. Possible values
#'   include "gaussian" (linear regression), "binomial" (logistic regression),
#'   "multinomial" (multinomial regression), and "cox" (survival analysis). This
#'   argument is only used with implementation="glmnet", or with functions using
#'   the family argument in the same way (see example below).
#' @param implementation name of the function to use for variable selection.
#'   With implementation="glmnet", the function \code{\link[glmnet]{glmnet}} is called.
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
#' @param n_cores number of cores to use for parallel computing. Only available
#'   on Unix systems.
#' @param verbose logical indicating if a message with minimum and maximum
#'   numbers of selected variables on one instance of resampled data should be
#'   printed.
#' @param ... additional parameters passed to the functions provided in
#'   "implementation" or "resampling".
#'
#' @return A list with: \item{S}{a matrix of the best stability scores for
#'   different penalty parameters. Rows correspond to different sets of penalty
#'   parameters, (values are stored in the output "Lambda").} \item{Lambda}{a
#'   matrix with different penalty parameters as rows.} \item{Q}{a matrix of
#'   average numbers of features selected by the underlying algorihm for
#'   different penalty parameters (as rows).} \item{Q_s}{a matrix of calibrated
#'   numbers of stable features for different penalty parameters (as rows).}
#'   \item{P}{a matrix of calibrated thresholds in selection proportions for
#'   different penalty parameters (as rows).} \item{PFER}{a matrix of computed
#'   upper-bounds in PFER of calibrated stability selection models for different
#'   penalty parameters.} \item{PFER}{a matrix of computed upper-bounds in FDP
#'   of calibrated stability selection models for different penalty parameters.}
#'   \item{S_2d}{a matrix of stability scores obtained with different
#'   combinations of parameters. Rows correspond to different penalty parameters
#'   and columns correspond to different tresholds in selection proportions.}
#'   \item{selprop}{a matrix of selection proportions. Columns correspond to
#'   predictors. Rows correspond to different penalty parameters.}
#'   \item{Beta}{an array of model coefficients. Rows correspond to different
#'   model parameters. Columns correspond to predictors. Indices along the third
#'   dimension correspond to different resampling iterations. With multivariate
#'   outcomes, indices along the fourth dimension correspond to outcome-specific
#'   coefficients.} \item{method}{a list with input values for the arguments
#'   "implementation", "family", "resampling" and "PFER_method".} \item{param}{a
#'   list with input values for the arguments "K", "pi_list", "tau", "n_cat",
#'   "pk", "PFER_thr", "FDP_thr", "seed", "xdata" and "ydata".}
#'
#' @seealso \code{\link{LambdaGridRegression}}, \code{\link{Resample}}
#'
#' @examples
#' # Variable selection in linear regression
#' simul=SimulateRegression(n=100, pk=20, family="gaussian") # simulated data
#' out=VariableSelection(xdata=simul$X, ydata=simul$Y, family="gaussian")
#'
#' @export
VariableSelection=function(xdata, ydata=NULL, Lambda=NULL, pi_list=seq(0.6,0.9,by=0.01),
                           K=100, tau=0.5, seed=1, n_cat=3,
                           family="gaussian", implementation="glmnet",
                           resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                           Lambda_cardinal=100,
                           n_cores=1, verbose=TRUE, ...){
  # Object preparation, error and warning messages
  CheckInputRegression(xdata=xdata, ydata=ydata, Lambda=Lambda, pi_list=pi_list,
                       K=K, tau=tau, seed=seed, n_cat=n_cat,
                       family=family, implementation=implementation,
                       resampling=resampling, PFER_method=PFER_method,
                       PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                       Lambda_cardinal=Lambda_cardinal,
                       verbose=verbose)

  if (is.null(Lambda)){
    # Defining grid of lambda values (using glmnet implementation)
    Lambda=LambdaGridRegression(xdata=xdata, ydata=ydata, tau=tau, seed=seed,
                                family=family, implementation=implementation,
                                resampling=resampling,
                                Lambda_cardinal=Lambda_cardinal, verbose=FALSE, ...)
  }

  # Stability selection and score
  mypar=parallel::mclapply(X=1:n_cores, FUN=function(k){
    return(SerialRegression(xdata=xdata, ydata=ydata, Lambda=Lambda, pi_list=pi_list,
                                    K=ceiling(K/n_cores), tau=tau, seed=as.numeric(paste0(seed, k)), n_cat=n_cat,
                                    family=family, implementation=implementation, resampling=resampling,
                                    PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                                    verbose=verbose, ...))
  })

  # Combining the outputs from parallel iterations
  out=mypar[[1]]
  if (n_cores>1){
    for (i in 2:length(mypar)){
      out=do.call(Combine, list(stability1=out, stability2=mypar[[2]], graph=FALSE))
    }
  }

  return(out)
}


#' Stability selection in regression (internal)
#'
#' This function can be used to estimate a
#' stability selection regression model
#' using a serial implementation and
#' when the grid of penalty parameters is provided
#' (for internal use only).
#'
#' @param xdata matrix of predictors with observations as rows and variables as columns.
#' @param ydata vector or matrix of outcome(s).
#' @param Lambda matrix of parameters controlling the underlying
#' feature selection algorithm specified in "implementation".
#' With implementation="glmnet", these are penalty parameters
#' controlling the regularised model.
#' @param pi_list grid of values for the threshold in selection proportion.
#' With n_cat=3, these values must be between 0.5 and 1.
#' With n_cat=2, these values must be between 0 and 1.
#' @param K number of resampling iterations.
#' @param tau subsample size. Only used with resampling="subsampling".
#' @param seed value of the seed to use to ensure reproducibility.
#' @param n_cat number of categories used to compute the stability score.
#' Possible values are 2 or 3.
#' @param family type of regression model. This argument is defined as in the
#' \code{\link[glmnet]{glmnet}} function from the glmnet package. Possible values include
#' "gaussian" (linear regression), "binomial" (logistic regression),
#' "multinomial" (multinomial regression), and "cox" (survival analysis).
#' This argument is only used with implementation="glmnet", or with functions
#' using the family argument in the same way (see example below).
#' @param implementation name of the function to use for variable selection.
#' With implementation="glmnet", the function \code{\link[glmnet]{glmnet}}
#' is called. Alternatively, this argument can be
#' a character string indicating the name of a function.
#' The function provided must use arguments called "x", "y", "lambda" and "family"
#' and return matrices of model coefficients (see \code{\link{SelectionAlgo}}).
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
#' @param FDP_thr threshold in the expected proportion of falsely selected features
#' (or False Discovery Proportion, FDP)
#' for constrained calibration by error control.
#' With PFER_thr=Inf and FDP_thr=Inf, unconstrained calibration is used.
#' @param verbose logical indicating if a message with minimum and maximum
#' numbers of selected variables on one instance of resampled data should be printed.
#' @param ... additional parameters passed to the functions provided in
#' "implementation" or "resampling".
#'
#' @return A list with:
#' \item{S}{a matrix of
#' the best stability scores
#' for different penalty parameters.
#' Rows correspond to different sets of penalty parameters,
#' (values are stored in the output "Lambda").}
#' \item{Lambda}{a matrix with different penalty parameters as rows.}
#' \item{Q}{a matrix of
#' average numbers of features
#' selected by the underlying algorihm
#' for different penalty parameters (as rows).}
#' \item{Q_s}{a matrix of
#' calibrated numbers of stable features
#' for different penalty parameters (as rows).}
#' \item{P}{a matrix of
#' calibrated thresholds in selection proportions
#' for different penalty parameters (as rows).}
#' \item{PFER}{a matrix of
#' computed upper-bounds in PFER of
#' calibrated stability selection models
#' for different penalty parameters.}
#' \item{PFER}{a matrix of
#' computed upper-bounds in FDP of
#' calibrated stability selection models
#' for different penalty parameters.}
#' \item{S_2d}{a matrix of
#' stability scores obtained
#' with different combinations of parameters.
#' Rows correspond to different penalty parameters and
#' columns correspond to different tresholds in selection proportions.}
#' \item{selprop}{a matrix of selection proportions.
#' Columns correspond to predictors.
#' Rows correspond to
#' different penalty parameters.}
#' \item{Beta}{an array of model coefficients.
#' Rows correspond to
#' different model parameters.
#' Columns correspond to predictors.
#' Indices along the third dimension correspond to
#' different resampling iterations.}
#' \item{method}{a list with input values for the arguments
#' "implementation", "family", "resampling" and "PFER_method".}
#' \item{param}{a list with input values for the arguments
#' "K", "pi_list", "tau", "n_cat", "pk", "PFER_thr", "FDP_thr",
#' "seed", "xdata" and "ydata".}
#'
#' @keywords internal
SerialRegression=function(xdata, ydata=NULL, Lambda, pi_list=seq(0.6,0.9,by=0.01),
                                  K=100, tau=0.5, seed=1, n_cat=3,
                                  family="gaussian", implementation="glmnet",
                                  resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                                  verbose=TRUE, ...){
  # Defining K if using complementary pairs (SS)
  if (PFER_method=="SS"){
    K=ceiling(K/2)*2
    tau=0.5
  }

  # Initialising objects to be filled
  N=N_block=ncol(xdata)
  Beta=array(0, dim=c(nrow(Lambda), ncol(xdata), K))
  rownames(Beta)=rownames(Lambda)
  colnames(Beta)=colnames(xdata)

  # Initialising the array with all beta coefficients
  s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
  Xsub = xdata[s,]
  Ysub = ydata[s,]
  mybeta=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
  if (length(dim(mybeta$beta_full))==2){
    Beta_full=array(0, dim=c(nrow(Lambda), dim(mybeta$beta_full)[2], K),
                    dimnames=list(rownames(Lambda), dimnames(mybeta$beta_full)[[2]], NULL))
  } else {
    if (length(dim(mybeta$beta_full))==3){
      Beta_full=array(0, dim=c(nrow(Lambda), dim(mybeta$beta_full)[2], K, dim(mybeta$beta_full)[3]),
                      dimnames=list(rownames(Lambda), dimnames(mybeta$beta_full)[[2]], NULL, dimnames(mybeta$beta_full)[[3]]))
    } else {
      stop(paste0("Invalid output from the variable selection function: ", implementation, "(). The output 'beta_full' must be an array with 2 or 3 dimensions."))
    }
  }

  # Computation of the selection proportions over Lambda
  if (verbose){
    pb=utils::txtProgressBar(style=3)
  }
  if (PFER_method=="MB"){
    for (k in 1:K){
      set.seed(k)
      s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      mybeta=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Resampling if model failed to converge
      while (is.infinite(mybeta$selected[1])){
        s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
        Xsub = xdata[s,]
        Ysub = ydata[s,]
        mybeta=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      }

      # Storing (one set of) beta coefficients, used to define set of selected variables
      Beta[rownames(mybeta$selected),colnames(mybeta$selected),k]=mybeta$selected

      # Storing all beta coefficients
      if (length(dim(Beta_full))==3){
        Beta_full[rownames(mybeta$beta_full),colnames(mybeta$beta_full),k]=mybeta$beta_full
      } else {
        Beta_full[rownames(mybeta$beta_full),colnames(mybeta$beta_full),k,]=mybeta$beta_full
      }

      if (verbose){
        utils::setTxtProgressBar(pb, k/K)
      }
    }

    # Computing the selection proportions
    bigstab=apply(Beta,c(1,2),FUN=function(x){sum(x!=0)})/K
  }

  if (PFER_method=="SS"){
    for (k in 1:ceiling(K/2)){
      set.seed(k)
      s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

      # First subset
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      mybeta1=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Complementary subset
      Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      mybeta2=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Resampling if model failed to converge
      while (is.infinite(mybeta1$selected[1])|is.infinite(mybeta2$selected[1])){
        s=Resample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

        # First subset
        Xsub = xdata[s,]
        Ysub = ydata[s,]
        mybeta=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

        # Complementary subset
        Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
        Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
        mybeta=SelectionAlgo(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      }

      # Storing beta coefficients from first set
      Beta[rownames(mybeta1$selected),colnames(mybeta1$selected),k]=mybeta1$selected

      # Storing all beta coefficients from first set
      if (length(dim(Beta_full))==3){
        Beta_full[rownames(mybeta1$beta_full),colnames(mybeta1$beta_full),k]=mybeta1$beta_full
      } else {
        Beta_full[rownames(mybeta1$beta_full),colnames(mybeta1$beta_full),k,]=mybeta1$beta_full
      }

      # Storing beta coefficients from complementary set
      Beta[rownames(mybeta2$selected),colnames(mybeta2$selected),ceiling(K/2)+k]=mybeta2$selected

      # Storing all beta coefficients from complementary set
      if (length(dim(Beta_full))==3){
        Beta_full[rownames(mybeta2$beta_full),colnames(mybeta2$beta_full),ceiling(K/2)+k]=mybeta2$beta_full
      } else {
        Beta_full[rownames(mybeta2$beta_full),colnames(mybeta2$beta_full),ceiling(K/2)+k,]=mybeta2$beta_full
      }

      if (verbose){
        utils::setTxtProgressBar(pb, 2*k/K)
      }
    }

    # Computing the simultaneous selection proportions
    bigstab=matrix(0, nrow=nrow(Beta), ncol=ncol(Beta))
    for (k in 1:ceiling(K/2)){
      A1=ifelse(Beta[,,k]!=0, yes=1, no=0)
      A2=ifelse(Beta[,,ceiling(K/2)+k]!=0, yes=1, no=0)
      A=A1+A2
      A=ifelse(A==2, yes=1, no=0)
      bigstab=bigstab+A
    }
    bigstab=bigstab/ceiling(K/2)
  }
  cat("\n")

  # Computation of the stability score over Lambda and pi_list
  if (K>2){
    metrics=StabilityMetrics(bigstab=bigstab, pk=NULL, pi_list=pi_list, K=K, n_cat=n_cat,
                           Sequential_template=NULL, graph=FALSE,
                           PFER_method=PFER_method, PFER_thr_blocks=PFER_thr, FDP_thr_blocks=FDP_thr)
    if (verbose){
      utils::setTxtProgressBar(pb, 1)
      cat("\n")
    }
  } else {
    Q=matrix(NA, nrow=nrow(Lambda), ncol=1)
    for (k in 1:nrow(Lambda)){
      q_block=sum(Beta[k,,1]!=0)
      Q[k,1]=round(q_block)
    }
  }
  Beta=Beta_full

  # Preparing outputs
  if (K>2){
    return(list(S=metrics$S, Lambda=Lambda,
                Q=metrics$Q, Q_s=metrics$Q_s, P=metrics$P,
                PFER=metrics$PFER, FDP=metrics$FDP,
                S_2d=metrics$S_2d, PFER_2d=metrics$PFER_2d, FDP_2d=metrics$FDP_2d,
                selprop=bigstab, Beta=Beta,
                methods=list(implementation=implementation, family=family, resampling=resampling, PFER_method=PFER_method),
                params=list(K=K, pi_list=pi_list, tau=tau, n_cat=n_cat, pk=ncol(xdata), PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed, xdata=xdata, ydata=ydata)))
  } else {
    return(list(Q=Q, Beta=Beta))
  }
}



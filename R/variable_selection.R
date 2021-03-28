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
                                resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                                Lambda_cardinal=Lambda_cardinal, verbose=FALSE, ...)
  }

  # Stability selection and score
  mypar=parallel::mclapply(X=1:n_cores, FUN=function(k){
    return(StabilityCalibRegression(xdata=xdata, ydata=ydata, Lambda=Lambda, pi_list=pi_list,
                                    K=ceiling(K/n_cores), tau=tau, seed=k, n_cat=n_cat,
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


# Not to be accessible
StabilityCalibRegression=function(xdata, ydata=NULL, Lambda, pi_list=seq(0.6,0.9,by=0.01),
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
  s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
  Xsub = xdata[s,]
  Ysub = ydata[s,]
  mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
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
      s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Resampling if model failed to converge
      while (is.infinite(mybeta$beta[1])){
        s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)
        Xsub = xdata[s,]
        Ysub = ydata[s,]
        mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      }

      # Storing (one set of) beta coefficients, used to define set of selected variables
      Beta[rownames(mybeta$beta),colnames(mybeta$beta),k]=mybeta$beta

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
      s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

      # First subset
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      mybeta1=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Complementary subset
      Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      mybeta2=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

      # Resampling if model failed to converge
      while (is.infinite(mybeta1$beta[1])|is.infinite(mybeta2$beta[1])){
        s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

        # First subset
        Xsub = xdata[s,]
        Ysub = ydata[s,]
        mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)

        # Complementary subset
        Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
        Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
        mybeta=SelectionFunction(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, implementation=implementation, ...)
      }

      # Storing beta coefficients from first set
      Beta[rownames(mybeta1$beta),colnames(mybeta1$beta),k]=mybeta1$beta

      # Storing all beta coefficients from first set
      if (length(dim(Beta_full))==3){
        Beta_full[rownames(mybeta1$beta_full),colnames(mybeta1$beta_full),k]=mybeta1$beta_full
      } else {
        Beta_full[rownames(mybeta1$beta_full),colnames(mybeta1$beta_full),k,]=mybeta1$beta_full
      }

      # Storing beta coefficients from complementary set
      Beta[rownames(mybeta2$beta),colnames(mybeta2$beta),ceiling(K/2)+k]=mybeta2$beta

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
    metrics=ComputeMetrics(bigstab=bigstab, pk=NULL, pi_list=pi_list, K=K, n_cat=n_cat,
                           Sequential_template=NULL, graph=FALSE,
                           PFER_method=PFER_method, PFER_thr_blocks=PFER_thr, FDP_thr_blocks=FDP_thr)
    if (verbose){
      utils::setTxtProgressBar(pb, 1)
      cat("\n")
    }
  } else {
    Q=matrix(NA, nrow=nrow(Lambda), ncol=1)
    PFER=matrix(NA, nrow=nrow(Lambda), ncol=length(pi_list))
    for (k in 1:nrow(Lambda)){
      q_block=sum(Beta[k,,1]!=0)
      Q[k,1]=round(q_block)
      for (j in 1:length(pi_list)){
        pi=pi_list[j]
        PFER[k,j]=ComputePFER(q=q_block,pi=pi,N=N_block)
      }
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
    return(list(Q=Q, Beta=Beta, PFER_2d=PFER))
  }
}



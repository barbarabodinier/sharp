LambdaGridRegression=function(xdata, ydata, tau=0.5, seed=1,
                              family="gaussian", implementation="glmnet",
                              resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                              Lambda_cardinal=100,
                              verbose=TRUE, ...){
  # Object preparation, error and warning messages
  Lambda=NULL
  pi_list=seq(0.6,0.9,by=0.01)
  K=100
  n_cat=3
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
  s=GetSubsample(data=ydata, family=family, tau=tau, resampling=resampling, ...)

  # Getting upperbound of Lambda
  set.seed(1) # To keep to allow for reproducible parallelisation
  if (implementation=="glmnet"){
    mycv=glmnet::cv.glmnet(x=xdata[s,], y=ydata[s,], family=family, ...)
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
    if (implementation=="gglasso"){
      print(paste("Minimum number of selected groups:", min(mycv$nzero)))
      print(paste("Maximum number of selected groups:", max(mycv$nzero)))
    }
  }

  # Creating a grid of lambda values from min and max
  Lambda=cbind(GetLambdaPath(lmax=max(mycv$lambda), lmin=min(mycv$lambda), cardinal=Lambda_cardinal))
  Lambda=as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))

  return(Lambda)
}


#' @export
LambdaGridNetwork=function(data, pk=NULL, lambda_other_blocks=NULL, K=100, tau=0.5, n_cat=3,
                           implementation="glassoFast", start="cold", scale=TRUE,
                           resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                           Lambda_cardinal=50, lambda_max=NULL, lambda_path_factor=0.001, max_density=0.5,
                           verbose=TRUE, ...){
  # K to keep for PFER computations with PFER_method set to "SS"
  # Error and warning messages
  bigblocks=bigblocks_vect=blocks=N_blocks=nblocks=PFER_thr_blocks=FDP_thr_blocks=NULL
  Lambda=NULL
  seed=1 # To keep to allow for reproducible parallelisation
  pi_list=seq(0.6,0.9,by=0.01)
  CheckInputNetwork(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=lambda_other_blocks,
                    pi_list=pi_list, K=K, tau=tau, seed=seed, n_cat=n_cat,
                    implementation=implementation, start=start, scale=scale,
                    resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                    Lambda_cardinal=Lambda_cardinal,
                    lambda_max=lambda_max, lambda_path_factor=lambda_path_factor, max_density=max_density,
                    verbose=verbose)
  rm(Lambda)
  pi_list=0.75 # only used for screening

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
        Lambda=cbind(Lambda, GetLambdaPath(lambda_max, lmin[b], cardinal=Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[l,,drop=FALSE]
        myscreen=StabilityCalibNetwork(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1,
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
        Lambda=cbind(Lambda, GetLambdaPath(lambda_max, lmin[b], cardinal=Lambda_cardinal))
      }

      # Initialisation of the smallest lambda
      lmin=Lambda[2,]
      l=1
      while (l<nrow(Lambda)){
        if (is.null(lambda_other_blocks)){
          ldense=lmin
        }
        tmpLambda=Lambda[l,,drop=FALSE]
        myscreen=StabilityCalibNetwork(data=data, pk=pk, Lambda=tmpLambda, lambda_other_blocks=ldense, pi_list=pi_list, K=1,
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
            mytmplist=c(mytmplist, ComputePFER(q=myscreen$Q[b,b], pi=pi, N=N_blocks[b], K=K, PFER_method=PFER_method))
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
    Lambda=cbind(Lambda, GetLambdaPath(lambda_max, lmin[b], cardinal=Lambda_cardinal))
  }
  Lambda=as.matrix(stats::na.exclude(Lambda))
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))

  return(Lambda)
}


GetLambdaPath=function(lmax, lmin, cardinal=100){
  return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}


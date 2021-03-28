GraphicalModel=function(data, pk=NULL, Lambda=NULL, lambda_other_blocks=NULL,
                        pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, n_cat=3,
                        implementation="glassoFast", start="warm", scale=TRUE,
                        resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                        Lambda_cardinal=50, lambda_max=NULL, lambda_path_factor=0.001, max_density=0.5,
                        n_cores=1, verbose=TRUE, ...){
  # Definition of the type of approach (single or multi-block)
  if (is.null(pk)){
    pk=ncol(data)
  }
  if (length(pk)>1){
    calibration="multi-block"
  } else {
    calibration="single-block"
  }
  if (verbose){
    print(paste("Starting", calibration, "calibration..."))
  }

  # Error and warning messages
  bigblocks=bigblocks_vect=blocks=N_blocks=nblocks=PFER_thr_blocks=FDP_thr_blocks=NULL
  CheckInputNetwork(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=lambda_other_blocks,
                    pi_list=pi_list, K=K, tau=tau, seed=seed, n_cat=n_cat,
                    implementation=implementation, start=start, scale=scale,
                    resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                    Lambda_cardinal=Lambda_cardinal,
                    lambda_max=lambda_max, lambda_path_factor=lambda_path_factor, max_density=max_density,
                    verbose=verbose)

  # Launching stability selection and calibration
  if (is.null(Lambda)){
    # Defining a broad grid of lambda values
    if (verbose){
      print("Defining the grid of lambda values...")
    }
    Lambda=LambdaGridNetwork(data=data, pk=pk, lambda_other_blocks=lambda_other_blocks, tau=tau,
                             implementation=implementation, start="cold", scale=scale,
                             resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                             lambda_max=lambda_max, lambda_path_factor=lambda_path_factor, max_density=max_density,
                             Lambda_cardinal=Lambda_cardinal, verbose=FALSE, ...)
  }

  # Check if parallelisation is possible (forking)
  if (.Platform$OS.type!="unix"){
    if (n_cores>1){
      warning("Invalid input for argument 'n_cores'. Parallelisation relies on forking, it is not available on Windows.")
    }
    n_cores=1
  }

  # Stability selection and score
  mypar=parallel::mclapply(X=1:n_cores, FUN=function(k){
    return(StabilityCalibNetwork(data=data, pk=pk, Lambda=Lambda, lambda_other_blocks=lambda_other_blocks,
                                 pi_list=pi_list, K=ceiling(K/n_cores), tau=tau, seed=k, n_cat=n_cat,
                                 implementation=implementation, start=start, scale=scale,
                                 resampling=resampling, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr,
                                 verbose=verbose, ...))
  })

  # Combining the outputs from parallel iterations
  out=mypar[[1]]
  if (n_cores>1){
    for (i in 2:length(mypar)){
      out=do.call(Combine, list(stability1=out, stability2=mypar[[2]], graph=TRUE))
    }
  }

  if (verbose){
    cat("\n")
    print("Visited Q:")
    if (nrow(out$Q)>15){
      print(utils::head(out$Q))
      print("[...]")
      print(utils::tail(out$Q))
    } else {
      print(out$Q)
    }
  }

  return(out)
}


StabilityCalibNetwork=function(data, pk=NULL, Lambda, lambda_other_blocks=NULL,
                               pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, n_cat=n_cat,
                               implementation="glassoFast", start="cold", scale=TRUE,
                               resampling="subsampling", PFER_method="MB", PFER_thr=Inf, FDP_thr=Inf,
                               verbose=TRUE, ...){
  # Marginal correlation to get sign of the relationship
  mycor_for_sign=stats::cor(data)

  # Defining K if using complementary pairs (SS)
  if (PFER_method=="SS"){
    K=ceiling(K/2)*2
    tau=0.5
  }

  # Creating matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)

  # Preparing the PFER and FDP thresholds
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

  # Re-formatting Lambda
  if (is.vector(Lambda)){
    grid=MakeBlockLambdaGrid(Lambda=Lambda, lambda_other_blocks=lambda_other_blocks)
    Lambda=grid$Lambda
    Sequential_template=grid$Sequential_template
  } else {
    grid=MakeBlockLambdaGrid(Lambda=Lambda, lambda_other_blocks=lambda_other_blocks)
    Lambda=grid$Lambda
    Sequential_template=grid$Sequential_template
  }

  # Showing the grid of (block-specific) lambda values
  if (verbose){
    print("Grid of lambda values:")
    if (ncol(Lambda)==1){
      print(as.vector(Lambda))
    } else {
      print(Lambda)
    }
  }

  # Initialising array of selection proportions
  bigstab=array(0, dim=c(ncol(data), ncol(data), nrow(Lambda)),
                dimnames=list(colnames(data), colnames(data), NULL))

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


  # Initialisation of the run
  if (verbose){
    pb=utils::txtProgressBar(style=3)
  }
  set.seed(seed)

  # Using MB formula of the PFER
  if (PFER_method=="MB"){
    for (i in 1:K){
      # Subsampling of the data
      s=GetSubsample(data=data, family=NULL, tau=tau, resampling=resampling, ...)
      data_sub=data[s,]

      # Estimation of the networks for different penalties
      A=NetworkFunction(x=data_sub, pk=pk, Lambda=Lambda, Sequential_template=Sequential_template,
                        scale=scale, implementation=implementation, start=start, ...)

      # Computing the selection counts
      for (k in 1:dim(A)[3]){
        bigstab[,,k]=bigstab[,,k]+A[,,k]
      }

      if (verbose){
        utils::setTxtProgressBar(pb, i/K)
      }
    }
    # Getting selection proportions from selection counts
    for (k in 1:dim(bigstab)[3]){
      bigstab[,,k]=bigstab[,,k]/K
      diag(bigstab[,,k])=0
    }
  }

  # Using complementary pairs and SS formula of the PFER
  if (PFER_method=="SS"){
    for (i in 1:ceiling(K/2)){
      # Sample 1
      s=GetSubsample(data=data, family=NULL, tau=tau, resampling=resampling, ...)
      data_sub=data[s,]

      # Estimation of the networks for different penalties
      A1=NetworkFunction(x=data_sub, pk=pk, Lambda=Lambda, Sequential_template=Sequential_template,
                        scale=scale, implementation=implementation, start=start, ...)

      # # Computing the selection counts
      # for (k in 1:dim(A)[3]){
      #   bigstab[,,k]=bigstab[,,k]+A[,,k]
      # }

      # Sample 2: everything not in sample 1
      data_sub=data[-s,]

      # Estimation of the networks for different penalties
      A2=NetworkFunction(x=data_sub, pk=pk, Lambda=Lambda, Sequential_template=Sequential_template,
                        scale=scale, implementation=implementation, start=start, ...)

      # Computing the simultaneous selection counts
      for (k in 1:dim(A1)[3]){
        A=ifelse((A1[,,k]+A2[,,k])==2, yes=1, no=0)
        bigstab[,,k]=bigstab[,,k]+A
      }

      if (verbose){
        utils::setTxtProgressBar(pb, i/ceiling(K/2))
      }
    }
    # Getting selection proportions from selection counts
    for (k in 1:dim(bigstab)[3]){
      bigstab[,,k]=bigstab[,,k]/ceiling(K/2)
      diag(bigstab[,,k])=0
    }
  }

  # Computation of the stability score
  if (K>2){
    metrics=ComputeMetrics(bigstab=bigstab, pk=pk, pi_list=pi_list, K=K, n_cat=n_cat,
                           Sequential_template=Sequential_template, graph=TRUE,
                           PFER_method=PFER_method, PFER_thr_blocks=PFER_thr_blocks, FDP_thr_blocks=FDP_thr_blocks)
    if (verbose){
      utils::setTxtProgressBar(pb, 1)
      cat("\n")
    }
  } else {
    # Initialising objects to be filled
    Q=matrix(NA,nrow=nrow(Lambda),ncol=nblocks)
    PFER=matrix(NA,nrow(Lambda),length(pi_list))
    for (k in 1:nrow(Lambda)){
      # Extracting corresponding selection proportions
      stab_iter=bigstab[,,k]

      # Getting number of selected variables per block
      for (block_id in 1:nblocks){
        stab_iter_block=stab_iter[(bigblocks==block_id)&(upper.tri(bigblocks))] # selection proportions in the block
        q_block=round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
        Q[k,block_id]=q_block
      }
    }
  }

  # Preparing outputs
  if (K>2){
    if (nblocks==1){
      return(list(S=metrics$S, Lambda=Lambda,
                  Q=metrics$Q, Q_s=metrics$Q_s, P=metrics$P,
                  PFER=metrics$PFER, FDP=metrics$FDP,
                  S_2d=metrics$S_2d, PFER_2d=metrics$PFER_2d, FDP_2d=metrics$FDP_2d,
                  selprop=bigstab, sign=sign(mycor_for_sign),
                  methods=list(implementation=implementation, start=start, resampling=resampling, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, n_cat=n_cat, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed,
                              lambda_other_blocks=lambda_other_blocks, Sequential_template=Sequential_template, data=data)))
    } else {
      return(list(S=metrics$S, Lambda=Lambda,
                  Q=metrics$Q, Q_s=metrics$Q_s, P=metrics$P,
                  PFER=metrics$PFER, FDP=metrics$FDP,
                  S_2d=metrics$S_2d,
                  selprop=bigstab, sign=sign(mycor_for_sign),
                  methods=list(implementation=implementation, start=start, resampling=resampling, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, n_cat=n_cat, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed,
                              lambda_other_blocks=lambda_other_blocks, Sequential_template=Sequential_template, data=data)))
    }
  } else {
    return(list(Q=Q))
  }
}


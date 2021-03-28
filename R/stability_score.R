ComputeMetrics=function(bigstab, pk=NULL, pi_list=seq(0.6,0.9,by=0.01),
                        K=100, n_cat=3,
                        PFER_method="MB", PFER_thr_blocks=Inf, FDP_thr_blocks=Inf,
                        Sequential_template=NULL, graph=TRUE){
  if (graph){
    nlambda=dim(bigstab)[3]
  } else {
    nlambda=nrow(bigstab)
  }

  if (is.null(pk)){
    pk=ncol(bigstab)
  }

  if (is.null(Sequential_template)){
    Sequential_template=matrix(TRUE, nrow=nlambda, ncol=1)
  }

  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)

  # Initialising objects to be filled
  Q=Q_s=P=matrix(NA,nrow=nlambda,ncol=nblocks)
  best_loglik=best_PFER=best_FDP=matrix(NA, nrow=nlambda, ncol=nblocks)
  if (nblocks==1){
    loglik=PFER=FDP=matrix(NA, ncol=length(pi_list), nrow=nlambda)
  } else {
    loglik=array(NA, dim=c(nlambda, length(pi_list), nblocks))
  }

  # Computing the metrics for each value of lambda
  for (k in 1:nlambda){
    # Extracting corresponding selection proportions
    if (graph){
      stab_iter=bigstab[,,k]
    } else {
      stab_iter=bigstab[k,]
    }

    # Computing stability score with block-specific pi
    for (block_id in 1:nblocks){
      if (Sequential_template[k,block_id]){
        if (graph){
          stab_iter_block=stab_iter[(bigblocks==block_id)&(upper.tri(bigblocks))] # selection proportions in the block
        } else {
          stab_iter_block=stab_iter
        }
        q_block=round(sum(stab_iter_block)) # average number of edges selected by the original procedure in the block
        Q[k,block_id]=q_block
        N_block=length(stab_iter_block) # maximum number of edges in the block
        tmp_loglik=tmp_PFERs=tmp_FDPs=rep(NA, length(pi_list))

        # Computing error rates and stability score for different values of pi
        for (j in 1:length(pi_list)){
          pi=pi_list[j]
          tmp_PFERs[j]=ComputePFER(q=q_block, pi=pi, N=N_block, K=K, PFER_method=PFER_method)
          tmp_FDPs[j]=ComputeFDP(PFER=tmp_PFERs[j], pi=pi, selprop=stab_iter_block)
          if ((tmp_PFERs[j]<=PFER_thr_blocks[block_id])&(tmp_FDPs[j]<=FDP_thr_blocks[block_id])){
            tmp_loglik[j]=StabilityScore(stab_iter=stab_iter_block, q=q_block, N=N_block, pi=pi, K=K, n_cat=n_cat)
          }
        }

        # Storing stability score in a matrix if only one block
        if (nblocks==1){
          loglik[k,]=tmp_loglik
          PFER[k,]=tmp_PFERs
          FDP[k,]=tmp_FDPs
        } else {
          loglik[k,,block_id]=tmp_loglik
        }

        # Keeping best stability score and other parameters at the max
        if (any(!is.na(tmp_loglik))){
          tmp_loglik[is.na(tmp_loglik)]=0
          myid=which.max(tmp_loglik)
          tmp_loglik[which(tmp_loglik==0)]=NA
          best_loglik[k,block_id]=tmp_loglik[myid]
          P[k,block_id]=pi_list[myid]
          Q_s[k,block_id]=sum(stab_iter_block>=pi_list[myid])
          best_PFER[k,block_id]=tmp_PFERs[myid]
          best_FDP[k,block_id]=tmp_FDPs[myid]
        }
      }
    }
  }
  best_loglik_blocks=best_loglik
  best_loglik=matrix(apply(best_loglik,1,sum), ncol=1)

  if (nblocks==1){
    return(list(S=best_loglik_blocks,
                Q=Q, Q_s=Q_s, P=P,
                PFER=best_PFER, FDP=best_FDP,
                S_2d=loglik, PFER_2d=PFER, FDP_2d=FDP))
  } else {
    return(list(S=best_loglik_blocks,
                Q=Q, Q_s=Q_s, P=P,
                PFER=best_PFER, FDP=best_FDP,
                S_2d=loglik))
  }
}


StabilityScore=function(stab_iter, q=NULL, N=NULL, pi, K, n_cat=3){
  # Preparing objects
  if (is.matrix(stab_iter)){
    stab_iter=stab_iter[upper.tri(stab_iter)]
  }

  # Computing the number of features (edges/variables)
  N=length(stab_iter)

  # Computing the average number of selected features
  if (is.null(q)){
    q=round(sum(stab_iter))
  }

  # Storing values of pi
  pi_list=pi

  # Loop over the values of pi
  score=rep(NA, length(pi_list))
  for (i in 1:length(pi_list)){
    pi=pi_list[i]

    # Computing the probabilities of being stable-in, stable-out or unstable under the null (uniform selection)
    p_vect=GetBinomialProbabilities(q, N, pi, K, n_cat=n_cat)

    # Computing the log-likelihood
    if (any(is.na(p_vect))){
      # Returning NA if not possible to compute (e.g. negative number of unstable features, as with pi<=0.5)
      l=NA
    } else {
      if (n_cat==2){
        S_0=sum(stab_iter<pi) # Number of not stably selected features
        S_1=sum(stab_iter>=pi) # Number of stably selected features

        # Checking consistency
        if (S_0+S_1!=N){
          stop(paste0("Inconsistency in number of edges \n S_0+S_1=", S_0+S_1, " instead of ", N))
        }

        # Log-likelihood
        l=S_0*p_vect$p_0+S_1*p_vect$p_1
      }

      if (n_cat==3){
        S_0=sum(stab_iter<=(1-pi)) # Number of stable-out features
        S_1=sum(stab_iter>=pi) # Number of stable-in features
        U=sum((stab_iter<pi)&(stab_iter>(1-pi))) # Number of unstable features

        # Checking consistency
        if (S_0+S_1+U!=N){
          stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0+S_1+U, " instead of ", N))
        }

        # Log-likelihood
        l=S_0*p_vect$p_1+U*p_vect$p_2+S_1*p_vect$p_3
      }

      # Re-formatting if infinite
      if (is.infinite(l)){
        l=NA
      }
    }

    # Getting the stability score
    score[i]=-l
  }

  return(score)
}


GetBinomialProbabilities=function(q, N, pi, K, n_cat=3){
  if (n_cat==2){
    # Definition of the threshold in selection counts
    thr=round(K*pi) # Threshold above (>=) which the feature is stably selected

    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_0=stats::pbinom(thr-1, size=K, prob=q/N, log.p=TRUE) # proportion < pi

    # Probability of observing a selection count above thr_up under the null
    p_1=stats::pbinom(thr-1, size=K, prob=q/N, lower.tail=FALSE, log.p=TRUE) # proportion >= pi

    # Checking consistency between the three computed probabilities (should sum to 1)
    if (abs(exp(p_0)+exp(p_1)-1)>1e-3){
      print(paste("N:", N))
      print(paste("q:", q))
      print(paste("K:", K))
      print(paste("pi:", pi))
      stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_0+p_1=", exp(p_0)+exp(p_1)))
    }

    # Output the two probabilities under the assumption of uniform selection procedure
    return(list(p_0=p_0, p_1=p_1))
  }

  if (n_cat==3){
    # Definition of the two thresholds in selection counts
    thr_down=round(K*(1-pi)) # Threshold below (<=) which the feature is stable-out
    thr_up=round(K*pi) # Threshold above (>=) which the feature is stable-in

    # Probability of observing a selection count below thr_down under the null (uniform selection)
    p_1=stats::pbinom(thr_down, size=K, prob=q/N, log.p=TRUE) # proportion <= (1-pi)

    # Probability of observing a selection count between thr_down and thr_up under the null
    if ((thr_down)>=(thr_up-1)){
      # Not possible to compute (i.e. negative number of unstable features)
      p_2=NA
    } else {
      # Using cumulative probabilities
      p_2=log(stats::pbinom(thr_up-1, size=K, prob=q/N)-stats::pbinom(thr_down, size=K, prob=q/N)) # 1-pi < proportion < pi

      # Using sum of probabilities (should not be necessary)
      if (is.infinite(p_2)|is.na(p_2)){
        p_2=0
        for (i in seq(thr_down+1, thr_up-1)){
          p_2=p_2+stats::dbinom(i, size=K, prob=q/N)
        }
        p_2=log(p_2)
      }
    }

    # Probability of observing a selection count above thr_up under the null
    p_3=stats::pbinom(thr_up-1, size=K, prob=q/N, lower.tail=FALSE, log.p=TRUE) # proportion >= pi

    # Checking consistency between the three computed probabilities (should sum to 1)
    if (!is.na(p_2)){
      if (abs(exp(p_1)+exp(p_2)+exp(p_3)-1)>1e-3){
        print(paste("N:", N))
        print(paste("q:", q))
        print(paste("K:", K))
        print(paste("pi:", pi))
        stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1)+exp(p_2)+exp(p_3)))
      }
    }

    # Output the three probabilities under the assumption of uniform selection procedure
    return(list(p_1=p_1, p_2=p_2, p_3=p_3))
  }
}


Combine=function(stability1, stability2, graph=TRUE){
  if (any(stability1$Lambda!=stability2$Lambda)){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different Lambdas.")
  }
  if (any(do.call(c,stability1$methods)!=do.call(c,stability2$methods))){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (any(do.call(c,stability1$params[c("pk","tau","PFER_thr","FDP_thr")])!=do.call(c,stability2$params[c("pk","tau","PFER_thr","FDP_thr")]))){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were constructed using different methods.")
  }
  if (stability1$params$seed==stability2$params$seed){
    warning("Arguments 'stability1' and 'stability2' were obtained using the same seed.")
  }
  if (any(stability1$params$data!=stability2$params$data)){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$params$xdata!=stability2$params$xdata)){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$params$ydata!=stability2$params$ydata)){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }
  if (any(stability1$sign!=stability2$sign)){
    stop("Arguments 'stability1' and 'stability2' are not compatible. They were obtained from different datasets.")
  }

  # Extracting the parameters
  Lambda=stability1$Lambda
  mysign=stability1$sign
  pk=stability1$params$pk
  pi_list=stability1$params$pi_list
  n_cat=stability1$params$n_cat
  Sequential_template=stability1$params$Sequential_template
  PFER_method=stability1$methods$PFER_method
  PFER_thr=stability1$params$PFER_thr
  FDP_thr=stability1$params$FDP_thr
  mymethods=stability1$methods
  myparams=stability1$params

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

  # Computing the total number of iterations
  K=stability1$params$K+stability2$params$K
  myparams$K=K
  myparams$seed=Inf

  # Computing selection propotions
  bigstab=array(NA, dim=dim(stability1$selprop), dimnames=dimnames(stability1$selprop))
  if (graph){
    for (k in 1:dim(stability1$selprop)[3]){
      bigstab[,,k]=(stability1$selprop[,,k]*stability1$params$K+stability2$selprop[,,k]*stability2$params$K)/K
    }
  } else {
    for (k in 1:nrow(stability1$selprop)){
      bigstab[k,]=(stability1$selprop[k,]*stability1$params$K+stability2$selprop[k,]*stability2$params$K)/K
    }
  }

  # Computation of the stability score
  metrics=ComputeMetrics(bigstab=bigstab, pk=pk, pi_list=pi_list, K=K, n_cat=n_cat,
                         Sequential_template=Sequential_template, graph=graph,
                         PFER_method=PFER_method, PFER_thr_blocks=PFER_thr_blocks, FDP_thr_blocks=FDP_thr_blocks)

  if (nblocks==1){
    return(list(S=metrics$S, Lambda=Lambda,
                Q=metrics$Q, Q_s=metrics$Q_s, P=metrics$P,
                PFER=metrics$PFER, FDP=metrics$FDP,
                S_2d=metrics$S_2d, PFER_2d=metrics$PFER_2d, FDP_2d=metrics$FDP_2d,
                selprop=bigstab, sign=mysign,
                methods=mymethods,
                params=myparams))
  } else {
    return(list(S=metrics$S, Lambda=Lambda,
                Q=metrics$Q, Q_s=metrics$Q_s, P=metrics$P,
                PFER=metrics$PFER, FDP=metrics$FDP,
                S_2d=metrics$S_2d,
                selprop=bigstab, sign=mysign,
                methods=mymethods,
                params=myparams))
  }
}


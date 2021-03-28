SelectionFunction=function(x, y, lambda, family, implementation="glmnet", ...){
  # Making sure none of the variables has a null standard deviation
  mysd=apply(x,2,stats::sd)
  if (any(mysd==0)){
    for (k in which(mysd==0)){
      x[,k]=x[,k]+stats::rnorm(n=nrow(x), sd=min(mysd[mysd!=0])/100)
    }
  }
  x=scale(x)

  if (implementation=="glmnet"){
    # Running the regression
    if (family=="multinomial"){
      mymodel=glmnet::glmnet(x=x, y=y, lambda=lambda, family=family, type.multinomial="grouped", ...)
    } else {
      mymodel=glmnet::glmnet(x=x, y=y, lambda=lambda, family=family, ...)
    }

    if (!is.infinite(mymodel$lambda[1])){
      # Extracting and formatting the beta coefficients
      if (!family%in%c("mgaussian", "multinomial")){
        mybeta=stats::coef(mymodel)
        mybeta=t(as.matrix(mybeta))
        mybeta=mybeta[,colnames(x)] # removing the intercept if included

        # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
        if (any(mysd==0)){
          mybeta[,which(mysd==0)]=0
        }

        # Preparing the outputs
        beta=beta_full=mybeta
      } else {
        mybeta=array(NA, dim=c(length(lambda), ncol(x), ncol(y)),
                     dimnames=list(paste0("s",0:(length(lambda)-1)), colnames(x), colnames(y)))
        for (y_id in 1:ncol(y)){
          tmpbeta=stats::coef(mymodel)[[y_id]]
          tmpbeta=t(as.matrix(tmpbeta))
          tmpbeta=tmpbeta[,colnames(x),drop=FALSE] # removing the intercept if included
          mybeta[rownames(tmpbeta),colnames(tmpbeta),y_id]=tmpbeta

          # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
          if (any(mysd==0)){
            mybeta[,which(mysd==0),y_id]=0
          }
        }

        # Preparing the outputs
        beta_full=mybeta
        beta=mybeta[,,1,drop=FALSE]
      }
    } else {
      # Returning infinite beta is the model failed
      beta=beta_full=Inf
    }
  } else {
    # Applying user-defined function for variable selection
    mybeta=do.call(get(implementation), args=list(x=x, y=y, lambda=lambda, family=family, ...))
    beta=mybeta$beta
    beta_full=mybeta$beta_full

    # Setting the beta coefficient to zero for predictors with always the same value (null standard deviation)
    if (any(mysd==0)){
      if (length(dim(beta_full))==2){
        beta[,which(mysd==0)]=0
        beta_full[,which(mysd==0)]=0
      }
      if (length(dim(beta_full))==3){
        beta[,which(mysd==0)]=0
        beta_full[,which(mysd==0),]=0
      }
    }
  }

  return(list(beta=beta, beta_full=beta_full))
}


NetworkFunction=function(x, pk=NULL, Lambda, Sequential_template, scale=TRUE, implementation="glassoFast", start="cold", ...){
  if (is.null(pk)){
    pk=ncol(x)
  }

  # Identifying potential variables with null standard deviation in the subsample
  mysd=apply(x,2,stats::sd)
  if (any(mysd==0)){
    for (k in which(mysd==0)){
      x[,k]=x[,k]+stats::rnorm(n=nrow(x), sd=min(mysd[mysd!=0])/100)
    }
  }

  # Create matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(as.vector(bigblocks_vect))
  names(N_blocks)=blocks
  nblocks=max(blocks)

  # Initialisation of array storing adjacency matrices
  adjacency=array(NA, dim=c(ncol(x), ncol(x), nrow(Lambda)))

  # Going over different (sets) of penalty parameters
  for (k in 1:nrow(Lambda)){
    # Creating penalisation matrix
    if (nblocks>1){
      lambdamat=bigblocks
      for (b in 1:nblocks){
        lambdamat[bigblocks==b]=Lambda[k,b]
      }
    } else {
      lambdamat=Lambda[k,1]
    }

    if (implementation=="glassoFast"){
      # Estimation of the covariance
      if (scale){
        cov_sub=stats::cor(x)
      } else {
        cov_sub=stats::cov(x)
      }

      # Estimation of the sparse inverse covariance
      if ((start=="warm")&(k!=1)){
        if (which(Sequential_template[k,])==which(Sequential_template[k-1,])){
          g_sub=glassoFast::glassoFast(S=cov_sub, rho=lambdamat,
                                       start="warm", w.init=sigma, wi.init=omega)
        } else {
          # Cold start if first iteration for the block
          g_sub=glassoFast::glassoFast(S=cov_sub, rho=lambdamat)
        }
      } else {
        g_sub=glassoFast::glassoFast(S=cov_sub, rho=lambdamat)
      }
      omega=g_sub$wi
      sigma=g_sub$w

      # Creating adjacency matrix
      A=ifelse(omega!=0, yes=1, no=0)
      A=A+t(A)
      A=ifelse(A!=0, yes=1, no=0)
    } else {
      A=do.call(get(implementation), args=list(x=x, lambda=lambdamat, scale=scale, ...))
    }

    # Ensuring that there is no edge for variables with always the same value (null standard deviation)
    if (any(mysd==0)){
        A[which(mysd==0),]=0
        A[,which(mysd==0)]=0
    }

    # Storing the estimated adjacency matrix
    adjacency[,,k]=A
  }

  return(adjacency)
}





SimulateGraph=function(n=100, pk=10, implementation="huge", topology="random", nu=0.1,
                       output_matrices=FALSE,
                       v_within=1, v_between=1,
                       pd_strategy="diagonally_dominant",
                       u=NULL, niter_max_u_grid=5, tolerance_u_grid=10, u_delta=5, ...){
  # pd_strategy: "diagonally_dominant" or "nonnegative_eigenvalues"
  # n_iter_max_u_grid: number of iterations where grid is refined
  # tolerance_u_grid: how far away from the extreme values
  # u_delta: difference in log10-scale between the smallest (or largest) u values in the grid used in the current iteration and the one built for next iteration

  # Defining grid of u values if not provided
  if (is.null(u)){
    u=10^-(seq(0,5,by=0.1))
    refining_u_grid=TRUE
    niter_max=5
    tolerance=10
  } else {
    refining_u_grid=FALSE
  }

  # Defining number of nodes
  p=sum(pk)

  # Creating matrix with block indices
  bigblocks=GetBlockMatrix(pk)
  bigblocks_vect=bigblocks[upper.tri(bigblocks)]
  N_blocks=unname(table(bigblocks_vect))
  block_ids=unique(as.vector(bigblocks))
  names(N_blocks)=block_ids
  nblocks=max(block_ids)

  # Building v matrix
  v_list=rep(NA,length(block_ids))
  v_list[diag(bigblocks)]=v_within
  v_list[is.na(v_list)]=v_between
  v=bigblocks
  for (k in block_ids){
    v[bigblocks==k]=v_list[k]
  }

  # Simulation of the adjacency matrix
  if (implementation=="huge"){
    theta=SimulateAdjacency(pk=pk, topology=topology, nu=nu, ...)
  } else {
    theta=do.call(get(implementation), args=list(pk=pk, topology=topology, nu=nu, ...))
  }

  # Ensuring that there is no self-loops
  diag(theta)=0

  # Setting variable names
  colnames(theta)=rownames(theta)=paste0("var",1:ncol(theta))

  # Filling off-diagonal entries of the precision matrix
  omega=theta*v

  # Calibrate u based on contrasts of the correlation matrix
  contrast=NULL
  for (u_value in u){
    omega_tmp=MakePositiveDefinite(omega=omega, u_value=u_value, pd_strategy=pd_strategy)
    C=stats::cov2cor(solve(omega_tmp))
    contrast=c(contrast, GetContrast(C))
  }

  # Avoiding extreme values in u grid if not provided by the user
  if (refining_u_grid){
    stop=0
    niter=1
    while (stop==0){
      niter=niter+1
      if (niter==niter_max_u_grid){
        stop=1
      }
      # Satisfied with calibrated u if the argmax is not too close to the boundaries (as defined from tolerance_u_grid)
      if (any(which(contrast==max(contrast))%in%seq(tolerance_u_grid,length(u)-tolerance_u_grid)==TRUE)){
        stop=1
      } else {
        # Adding smaller values of u
        if (any(which(contrast==max(contrast))%in%seq(1,tolerance_u_grid)==TRUE)){
          u=c(u, 10^-seq(min(-log10(u))-u_delta, min(-log10(u)), by=0.1))
        }

        # Adding larger values of u
        if (any(which(contrast==max(contrast))%in%seq(length(u)-tolerance_u_grid,length(u))==TRUE)){
          u=c(u, 10^-seq(max(-log10(u)), max(-log10(u)+u_delta), by=0.1))
        }

        # Sorting values in u
        u=sort(u, decreasing=TRUE)

        # Computing the contrast for all visited values of u
        contrast=NULL
        for (u_value in u){
          omega_tmp=MakePositiveDefinite(omega=omega, u_value=u_value, pd_strategy=pd_strategy)
          C=stats::cov2cor(solve(omega_tmp))
          contrast=c(contrast, GetContrast(C))
        }
      }
    }
  }

  # Computing calibrated precision matrix
  if (length(u)>1){
    u_value=u[which.max(contrast)]
    omega=MakePositiveDefinite(omega=omega, u_value=u_value, pd_strategy=pd_strategy)
  } else {
    omega=omega_tmp
  }

  # Simulating the sign of the correlations (uniform)
  sign_mat=matrix(0, nrow=nrow(omega), ncol=ncol(omega))
  sign_mat[upper.tri(sign_mat)]=sample(c(-1,1), size=sum(upper.tri(omega)), replace=TRUE)
  sign_mat=sign_mat+t(sign_mat)
  diag(sign_mat)=1
  omega=omega*sign_mat

  # Computing the correlation matrix
  C=stats::cov2cor(solve(omega)) # true correlation matrix - called sigma in huge

  # Computing the partial correlation matrix
  if (output_matrices){
    phi=-stats::cov2cor(omega)+2*diag(ncol(omega))
  }

  # Simulating data from multivariate normal distribution
  x=MASS::mvrnorm(n, rep(0, p), C)
  colnames(x)=paste0("var",1:ncol(x))
  rownames(x)=paste0("obs",1:nrow(x))

  if (output_matrices){
    return(list(data=x, theta=theta,
                omega=omega, phi=phi, C=C,
                u=u_value, u_grid=u, contrast_path=contrast))
  } else {
    return(list(data=x, theta=theta))
  }
}


SimulateXY=function(n=100, pk=10, X=NULL, nu_pred=0.2, beta_set=c(-1,1), continuous=FALSE,
                    sd_pred_error=1, family="gaussian"){
  # Simulation of the predictors
  if (is.null(X)) {
    p=sum(pk)
    X=NULL
    for (k in 1:p){
      X=cbind(X, stats::rnorm(n, mean=0, sd=1))
    }
    X=scale(X)
  } else {
    n=nrow(X)
    p=ncol(X)
  }

  # Setting column names for predictors
  if (is.null(colnames(X))){
    colnames(X)=paste0("var", 1:ncol(X))
  }

  # Setting row names for predictors
  if (is.null(rownames(X))){
    rownames(X)=paste0("obs", 1:nrow(X))
  }

  # Getting the binary vector of true predictors
  theta_pred=stats::rbinom(p, size=1, prob=nu_pred)
  names(theta_pred)=colnames(X)

  # Simulating a vector of betas
  if (continuous){
    beta=stats::runif(p, min=min(beta_set), max=max(beta_set))
  } else {
    beta=base::sample(beta_set, size=p, replace=TRUE)
  }
  beta=beta*theta_pred

  # Computing the predicted values of Y
  Y_pred=X%*%beta

  # Introducing some centered gaussian error
  Y=Y_pred+stats::rnorm(n, mean=0, sd=sd_pred_error)

  # Compute binary outcome for logistic regression
  if (family=="binomial"){
    proba=1/(1+exp(-Y)) # inverse logit
    Y_bin=stats::rbinom(n, size=1, prob=proba)
  }

  # Return the simulated X and Y
  if (family=="binomial"){
    out=list(X=X, Y=Y_bin, proba=proba, logit_proba=Y, logit_proba_pred=Y_pred, theta_pred=theta_pred, beta=beta)
  } else {
    out=list(X=X, Y=Y, Y_pred=Y_pred, theta_pred=theta_pred, beta=beta)
  }

  return(out)
}


SimulateAdjacency=function(pk=10, topology="random", nu=0.1, ...){
  # Simulating the adjacency matrix using huge
  theta=as.matrix(huge::huge.generator(n=2, d=sum(pk), prob=nu,
                                       graph=topology, verbose=FALSE, ...)$theta)

  # Re-organising the variables to avoid having centrality related to variable ID (e.g. for scale-free models)
  ids=sample(ncol(theta))
  theta=theta[ids,ids]

  return(theta)
}


MakePositiveDefinite=function(omega, u_value=0.1, pd_strategy="diagonally_dominant"){
  # Adding a small number (u) to the diagonal
  if (pd_strategy=="diagonally_dominant"){
    diag(omega)=apply(abs(omega),1,sum)+u_value
  }

  # Ensuring positive eigenvalues
  if (pd_strategy=="nonnegative_eigenvalues"){
    diag(omega)=abs(min(eigen(omega)$values))+u_value # used in huge
  }

  return(omega)
}


GetContrast=function(mat, digits=3){
  return(length(unique(round(as.vector(abs(mat)), digits=digits))))
}


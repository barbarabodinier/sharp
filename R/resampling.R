GetSubsample=function(data, family=NULL, tau=0.5, resampling="subsampling", ...){
  # Preparing the data
  if (is.vector(data)){
    data=matrix(data, ncol=1)
  }

  if (!resampling%in%c("subsampling","bootstrap")){
    s=do.call(get(resampling), args=list(data=data, tau=tau, ...))
  } else {
    # Using or not replacement in resampling
    replacement=ifelse(resampling=="subsampling", yes=FALSE, no=TRUE)

    # Definition of the size of sub/bootstrap sample
    if (replacement){
      tau=1
    }

    # Resampling procedure
    if (!is.null(family)){
      # Resampling for regression models
      if (family%in%c("gaussian","poisson","mgaussian")){
        s=sample(nrow(data), size=tau*nrow(data), replace=replacement)
      }
      if (family=="binomial"){
        s0=sample(which(data=="0"), size=tau*sum(data=="0"), replace=replacement)
        s1=sample(which(data=="1"), size=tau*sum(data=="1"), replace=replacement)
        s=c(s0,s1)
      }
      if (family=="multinomial"){
        s=NULL
        for (mycat in levels(factor(data))){
          scat=sample(which(data==mycat), size=tau*sum(data==mycat), replace=replacement)
          s=c(s,scat)
        }
      }
      if (family=="cox"){
        s0=sample(which(data[,2]=="0"), size=tau*sum(data[,2]=="0"), replace=replacement)
        s1=sample(which(data[,2]=="1"), size=tau*sum(data[,2]=="1"), replace=replacement)
        s=c(s0,s1)
      }
    } else {
      # Resampling for network models
      s=sample(1:nrow(data), size=tau*nrow(data), replace=replacement)
    }
  }
  return(s)
}



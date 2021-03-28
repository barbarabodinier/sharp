ComputePFER=function(q, pi, N, K, PFER_method="MB"){
  # Checking the inputs (PFER_method)
  PFER_method=as.character(PFER_method)
  if ((length(PFER_method)!=1)|(!PFER_method%in%c("MB","SS"))){
    stop("Invalid input for argument 'PFER_method'. Possible values are: 'MB' or 'SS'.")
  }

  if (pi>0.5){
    # Computing upper-bound of the PFER using approach proposed by MB
    if (PFER_method=="MB"){
      upperbound=1/(2*pi-1)*q^2/N
    }

    # Computing upper-bound of the PFER using approach proposed by SS
    if (PFER_method=="SS"){
      cutoff=pi
      B=ceiling(K/2)
      theta=q/N
      if (cutoff <= 3/4) {
        tmp=2*(2*cutoff-1-1/(2*B))
      } else {
        tmp=(1+1/B)/(4*(1-cutoff+1/(2*B)))
      }
      upperbound=q^2/N/tmp

      # Setting to Inf if "out of bounds"
      if ((cutoff<1/2+min(theta^2, 1/(2*B)+3/4*theta^2))|(cutoff>1)){
        upperbound=Inf
      }
    }
  } else {
    upperbound=Inf
  }

  # Re-formatting the upperbound
  if (is.na(upperbound)){
    upperbound=Inf
  }

  return(upperbound)
}


ComputeFDP=function(PFER, selprop, pi){
  # Formatting the selection proportions
  if (length(dim(selprop))==2){
    selprop=selprop[upper.tri(selprop)] # vector of selection proportions
  } else {
    selprop=selprop
  }

  # Computing the number of stable edges
  S=sum(selprop>=pi,na.rm=TRUE)

  # Computing the proportion of false discoveries among discoveries (False Discovery Proportion)
  if (S!=0){
    FDP=PFER/S
  } else {
    FDP=0
  }

  return(FDP)
}

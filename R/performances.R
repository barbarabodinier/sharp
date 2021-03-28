SelectionPerformance=function(theta, theta_star, pk=NULL, cor=NULL, thr=NULL){
  # Storing similarities/differences between estimated and true sets
  Asum=theta+2*theta_star

  # Extracting block-specific performances
  if (is.null(pk)){
    return(GetPerformance(Asum, cor=cor, thr=thr))
  } else {
    Asum_vect=Asum[upper.tri(Asum)]
    bigblocks=GetBlockMatrix(pk)
    bigblocks_vect=bigblocks[upper.tri(bigblocks)]
    if (!is.null(cor)){
      cor_vect=cor[upper.tri(cor)]
    } else {
      cor_vect=NULL
    }

    return(rbind(GetPerformance(Asum, cor=cor, thr=thr),
                 GetPerformance(Asum_vect[bigblocks_vect==1], cor=cor_vect[bigblocks_vect==1], thr=thr),
                 GetPerformance(Asum_vect[bigblocks_vect==2], cor=cor_vect[bigblocks_vect==2], thr=thr),
                 GetPerformance(Asum_vect[bigblocks_vect==3], cor=cor_vect[bigblocks_vect==3], thr=thr)))
  }
}


GetPerformance=function(Asum, cor=NULL, thr=0.5){
  # Asum is an adjacency matrix with 3 for TP, 2 for FN, 1 for FP, and 0 for TN

  # Preparing objects
  if (is.matrix(Asum)){
    p=ncol(Asum)
    N=p*(p-1)/2
    Asum=Asum[upper.tri(Asum)]
  } else {
    N=length(Asum)
  }

  # Computing the numbers of True/False Positives/Negatives
  TP=sum(Asum==3)
  FN=sum(Asum==2)
  FP=sum(Asum==1)
  TN=sum(Asum==0)

  # Separation being correlated and independent features based on a threshold in correlation
  if (!is.null(cor)){
    if (is.matrix(cor)){
      cor_vect=cor[upper.tri(cor)]
    } else {
      cor_vect=cor
    }
    FP_c=sum((Asum==1)&(abs(cor_vect)>=thr))
    FP_i=sum((Asum==1)&(abs(cor_vect)<thr))
  }

  # Computing performances in selection
  sensitivity=TP/(TP+FN)
  specificity=TN/(TN+FP)
  accuracy=(TP+TN)/N
  if (TP+FP>0){
    precision=TP/(TP+FP)
  } else {
    precision=0
  }
  if ((TP+FN)>0){
    recall=TP/(TP+FN)
  } else {
    recall=1
  }
  if ((precision>0)|(recall>0)){
    F1_score=2*precision*recall/(precision+recall)
  } else {
    F1_score=0
  }

  if (is.null(cor)){
    return(data.frame(TP=TP, FN=FN, FP=FP, TN=TN,
                      sensitivity=sensitivity, specificity=specificity,
                      accuracy=accuracy, precision=precision, recall=recall, F1_score=F1_score))
  } else {
    return(data.frame(TP=TP, FN=FN, FP=FP, TN=TN, FP_c=FP_c, FP_i=FP_i,
                      sensitivity=sensitivity, specificity=specificity,
                      accuracy=accuracy, precision=precision, recall=recall, F1_score=F1_score))
  }
}


SelectionPerformanceGraph=function(theta, theta_star, node_colour=NULL, plot=FALSE, filename=NULL, width=7, height=7,
                                   mycolours=c("grey95", "tomato", "forestgreen", "navy"),
                                   mylty=c(4,2,3,1), ...){
  # Storing similarities/differences between estimated and true sets
  Asum=theta+2*theta_star

  # Refining inputs
  names(mycolours)=names(mylty)=c("TN", "FP", "FN", "TP")
  if (is.null(node_colour)){
    node_colour=rep("black",ncol(Asum))
  }

  # Making consensus graph
  g=Graph(adjacency=ifelse(Asum!=0,yes=1,no=0), node_colour=node_colour, ...)

  # Formatting vertices
  igraph::V(g)$size=igraph::V(g)$size/3+1
  igraph::V(g)$label=rep("", length(igraph::V(g)$label))

  # Formatting edges
  myedgecolour=mycolours[Asum[igraph::get.edgelist(g)]+1]
  myedgelty=c(2, 3, 1)[Asum[igraph::get.edgelist(g)]]
  igraph::E(g)$color=myedgecolour
  igraph::E(g)$width=1
  igraph::E(g)$lty=myedgelty

  # Plotting graph
  if (plot|(!is.null(filename))){
    if (!is.null(filename)){
      grDevices::pdf(filename, width=width, height=height)
    }
    graphics::par(mar=rep(0,4))
    igraph::plot.igraph(g, layout=igraph::layout_with_kk(g))
    if (!is.null(filename)){
      grDevices::dev.off()
    }
  }

  # Returning output graph
  return(g)
}

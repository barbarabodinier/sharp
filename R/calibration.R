GetArgmaxId=function(stability=NULL, S=NULL){
  if ((is.null(stability))&(is.null(S))){
    stop("Invalid input. One of the two arguments has to be specified: 'stability' or 'S'.")
  }
  if (is.null(S)){
    argmax_id=matrix(NA, nrow=ncol(stability$Lambda), ncol=2)
    if (is.null(stability$params$lambda_other_blocks)&(length(stability$params$pk)>1)){
      id=which.max(apply(stability$S,1,sum,na.rm=TRUE))
      argmax_id[,1]=rep(id, nrow(argmax_id))
      for (block_id in 1:ncol(stability$Lambda)){
        if (!is.na(stability$P[id,block_id])){
          argmax_id[block_id,2]=which(stability$params$pi_list==stability$P[id,block_id])
        }
      }
    } else {
      for (block_id in 1:ncol(stability$Lambda)){
        if (ncol(stability$Lambda)==1){
          myS=stability$S
        } else {
          myS=stability$S[,block_id,drop=FALSE]
        }
        myS[is.na(myS)]=0
        myid=which.max(myS[,1])
        argmax_id[block_id,]=c(myid, which(stability$params$pi_list==stability$P[myid,block_id]))
      }
    }
  } else {
    argmax_id=matrix(NA, nrow=1, ncol=2)
    myS=apply(S,1,max,na.rm=TRUE)
    myS[is.na(myS)]=0
    myid=which.max(myS)
    argmax_id[1,]=c(myid, max(which(S[myid,]==myS[myid])))
  }
  colnames(argmax_id)=c("lambda_id","pi_id")
  return(argmax_id)
}


GetArgmax=function(stability){
  argmax=matrix(NA, nrow=ncol(stability$Lambda), ncol=2)
  if (is.null(stability$params$lambda_other_blocks)&(length(stability$params$pk)>1)){
    id=which.max(apply(stability$S,1,sum,na.rm=TRUE))
    argmax[,1]=stability$Lambda[id,]
    argmax[,2]=stability$P[id,]
  } else {
    for (block_id in 1:ncol(stability$Lambda)){
      if (ncol(stability$Lambda)==1){
        myS=stability$S
      } else {
        myS=stability$S[,block_id,drop=FALSE]
      }
      myS[is.na(myS)]=0
      myid=which.max(myS[,1])
      argmax[block_id,]=c(stability$Lambda[myid,block_id], stability$P[myid,block_id])
    }
  }
  colnames(argmax)=c("lambda","pi")
  return(argmax)
}


Adjacency=function(stability, argmax_id=NULL){
  A=matrix(0, ncol=ncol(stability$selprop), nrow=nrow(stability$selprop))
  bigblocks=GetBlockMatrix(stability$params$pk)
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(stability)
    argmax=GetArgmax(stability)
  } else {
    argmax=NULL
    for (block_id in 1:ncol(stability$Lambda)){
      argmax=rbind(argmax, c(stability$Lambda[argmax_id[block_id,1],block_id],
                             stability$params$pi_list[argmax_id[block_id,2]]))
    }
  }
  for (block_id in 1:ncol(stability$Lambda)){
    A_block=ifelse(stability$selprop[,,argmax_id[block_id,1]]>=argmax[block_id,2],1,0)
    A_block[lower.tri(A_block)]=0
    A_block=A_block+t(A_block) # for symmetry
    if (length(stability$params$pk)>1){
      A_block[bigblocks!=block_id]=0
    }
    A=A+A_block
  }
  A[is.na(A)]=0
  return(A)
}


SelectedVariables=function(stability, argmax_id=NULL){
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(stability)
  }
  stability_selected=ifelse(stability$selprop[argmax_id[1],]>=stability$params$pi_list[argmax_id[2]],1,0)
  return(stability_selected)
}


SelectionProportions=function(stability, argmax_id=NULL){
  if ("data"%in%names(stability$params)){
    out=SelectionProportionsNetwork(stability=stability, argmax_id=argmax_id)
  } else {
    out=SelectionProportionsRegression(stability=stability, argmax_id=argmax_id)
  }
  return(out)
}


SelectionProportionsNetwork=function(stability, argmax_id=NULL){
  A=matrix(0, ncol=ncol(stability$selprop), nrow=nrow(stability$selprop))
  bigblocks=GetBlockMatrix(stability$params$pk)
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(stability)
    argmax=GetArgmax(stability)
  } else {
    argmax=NULL
    for (block_id in 1:ncol(stability$Lambda)){
      argmax=rbind(argmax, c(stability$Lambda[argmax_id[block_id,1],],
                             stability$params$pi_list[argmax_id[block_id,2]]))
    }
  }
  for (block_id in 1:ncol(stability$Lambda)){
    A_block=stability$selprop[,,argmax_id[block_id,1]]
    A_block[lower.tri(A_block)]=0
    A_block=A_block+t(A_block) # for symmetry
    if (length(stability$params$pk)>1){
      A_block[bigblocks!=block_id]=0
    }
    A=A+A_block
  }
  return(A)
}


SelectionProportionsRegression=function(stability, argmax_id=NULL){
  if (is.null(argmax_id)){
    argmax_id=GetArgmaxId(stability)
    argmax=GetArgmax(stability)
  }
  m=stability$selprop[argmax_id[1],]
  calibrated_pi=stability$params$pi_list[argmax_id[2]]
  return(m)
}


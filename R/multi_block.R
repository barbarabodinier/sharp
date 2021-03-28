GetBlockMatrix=function(pk){
  nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
  blocks=matrix(NA, nrow=length(pk), ncol=length(pk))
  blocks[upper.tri(blocks, diag=TRUE)]=1:nblocks

  mybreaks=c(0,cumsum(pk))
  bigblocks=matrix(ncol=sum(pk), nrow=sum(pk))
  row_id_start=matrix(mybreaks[row(blocks)],ncol=length(pk))+1
  row_id_end=matrix(mybreaks[row(blocks)+1],ncol=length(pk))
  col_id_start=matrix(mybreaks[col(blocks)],ncol=length(pk))+1
  col_id_end=matrix(mybreaks[col(blocks)+1],ncol=length(pk))

  row_id_start=row_id_start[upper.tri(row_id_start, diag=TRUE)]
  row_id_end=row_id_end[upper.tri(row_id_end, diag=TRUE)]
  col_id_start=col_id_start[upper.tri(col_id_start, diag=TRUE)]
  col_id_end=col_id_end[upper.tri(col_id_end, diag=TRUE)]

  for (block_id in blocks[upper.tri(blocks, diag=TRUE)]){
    ids=rbind(expand.grid(row_id_start[block_id]:row_id_end[block_id],
                          col_id_start[block_id]:col_id_end[block_id]),
              expand.grid(col_id_start[block_id]:col_id_end[block_id],
                          row_id_start[block_id]:row_id_end[block_id]))
    bigblocks[as.matrix(ids)]=block_id
  }

  return(bigblocks)
}


MakeBlockLambdaGrid=function(Lambda, lambda_other_blocks=NULL){
  if ((is.null(lambda_other_blocks))&(!is.vector(Lambda))){
    Lambda_blocks=Lambda
    Sequential_template=matrix(TRUE, ncol=ncol(Lambda), nrow=nrow(Lambda))
  } else {
    # Create Lambda grid matrix with nblocks columns
    if (!is.null(lambda_other_blocks)){
      nblocks=length(lambda_other_blocks)
    } else {
      lambda_other_blocks=1
      nblocks=1
    }
    Lambda_blocks=NULL
    if (is.vector(Lambda)){
      Sequential_template=matrix(FALSE, nrow=nblocks*length(Lambda), ncol=nblocks)
    } else {
      Sequential_template=matrix(FALSE, nrow=nblocks*nrow(Lambda), ncol=nblocks)
    }
    for (block_id in 1:nblocks){
      if (!is.vector(Lambda)){
        tmpLambda=Lambda[,block_id]
      } else {
        tmpLambda=Lambda
      }
      Lambda_blocks=cbind(Lambda_blocks, rep(lambda_other_blocks[block_id], nblocks*length(tmpLambda)))
      Lambda_blocks[(length(tmpLambda)*(block_id-1)+1):(length(tmpLambda)*(block_id)), block_id]=tmpLambda
      Sequential_template[(length(tmpLambda)*(block_id-1)+1):(length(tmpLambda)*(block_id)), block_id]=TRUE
    }
  }

  return(list(Lambda=Lambda_blocks, Sequential_template=Sequential_template))
}



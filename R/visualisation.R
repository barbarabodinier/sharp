CalibrationPlot=function(stability, metric="both", block_id=NULL, filename=NULL, lines=TRUE,
                         xlab=expression(lambda), ylab=expression(pi), zlab=expression(italic(q)),
                         width=7, height=7, mar=NULL, mfrow=NULL, ...){
  # Extracting the number of blocks
  if (is.null(block_id)){
    bigblocks=GetBlockMatrix(stability$params$pk)
    bigblocks_vect=bigblocks[upper.tri(bigblocks)]
    N_blocks=unname(table(bigblocks_vect))
    blocks=unique(as.vector(bigblocks_vect))
    names(N_blocks)=blocks
    nblocks=max(blocks)
    block_id=1:nblocks
  } else {
    nblocks=1
  }

  # Saving as PDF
  if (!is.null(filename)){
    grDevices::pdf(filename, width=width, height=height)
  }

  if (metric=="both"){
    if (is.null(mar)){
      mar=c(7,5,7,7)
    }
    if (is.null(mfrow)){
      mfrow=c(1,nblocks)
    }
    graphics::par(mar=mar, mfrow=mfrow)

    for (b in block_id){
      # Extracting the stability scores
      if (length(stability$params$pk)==1){
        mat=stability$S_2d
        ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
        mat=mat[ids,,drop=FALSE]
      } else {
        mat=stability$S_2d[,,b]
        ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
        mat=mat[ids,,drop=FALSE]
      }
      mat=mat[,,drop=FALSE]
      colnames(mat)=stability$params$pi_list
      rownames(mat)=formatC(stability$Lambda[,b], format="e", digits=2)[ids]

      # Extracting corresponding numbers of selected variables (q)
      Q=stability$Q[,b]
      Q=Q[ids]

      # Heatmap representation
      Heatmap(mat[nrow(mat):1,ncol(mat):1])

      # Identifying best pair of parameters
      graphics::abline(h=which.min(abs(as.numeric(colnames(mat))-GetArgmax(stability)[b,2]))-0.5, lty=3)
      graphics::abline(v=nrow(mat)-which(stability$Lambda[ids,b]==GetArgmax(stability)[b,1])+0.5, lty=3)

      # Including axes
      graphics::axis(side=1, at=(1:nrow(mat))-0.5, las=2, labels=rev(rownames(mat)), ...)
      graphics::axis(side=2, at=(1:ncol(mat))-0.5, las=2,
                     labels=formatC(as.numeric(colnames(mat)), format="f", digits=2), ...)
      graphics::axis(side=3, at=(1:nrow(mat))-0.5, las=2,
                     labels=rev(formatC(Q, format="f", big.mark=",", digits=0)), ...)

      # Including axis labels
      graphics::mtext(text=xlab, side=1, line=5.2, cex=1.5)
      graphics::mtext(text=ylab, side=2, line=3.5, cex=1.5)
      graphics::mtext(text=zlab, side=3, line=3.5, cex=1.5)
    }
  } else {
    if (metric=="lambda"){
      if (is.null(mar)){
        mar=c(7,5,7,2)
      }
      if (is.null(mfrow)){
        mfrow=c(1,nblocks)
      }
      graphics::par(mar=mar, mfrow=mfrow)

      for (b in block_id){
        # Extracting the stability scores
        if (length(stability$params$pk)==1){
          mat=stability$S_2d
          ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
          mat=mat[ids,,drop=FALSE]
        } else {
          mat=stability$S_2d[,,b]
          ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
          mat=mat[ids,,drop=FALSE]
        }

        # Extracting the best stability score (with optimal pi) for each lambda value
        vect=apply(mat,1,max,na.rm=TRUE)

        # Extracting corresponding numbers of selected variables (q)
        Q=stability$Q[,b,drop=FALSE]
        Q=Q[ids]

        # Extracting corresponding lambda values
        Lambda=stability$Lambda[ids,b,drop=FALSE]

        # Re-ordering by decreasing lambda
        ids=sort.list(Lambda, decreasing=TRUE)
        Lambda=Lambda[ids]
        Q=Q[ids]
        vect=vect[ids]

        # Using input ylab if not as default
        ylab=ifelse(as.character(ylab)==as.character(expression(pi)), yes="Stability Score", no=ylab)

        # Making plot
        cex_points=0.7
        plot(Lambda, vect, pch=19, col="navy", cex=cex_points,
             xlab="", ylab=ylab, cex.lab=1.5, xaxt="n")
        graphics::abline(h=graphics::axTicks(side=2), lty=3, col="grey")
        graphics::abline(h=max(vect), lty=2, col="red")
        graphics::abline(v=which.max(vect), lty=2, col="red")
        graphics::points(Lambda, vect, pch=19, col="navy", cex=cex_points)
        if (lines){
          graphics::lines(Lambda, vect, col="navy")
        }

        # Adding x-axis and z-axis and their labels
        xseq=seq(1,length(Lambda),length.out=5)
        graphics::abline(v=Lambda[xseq], lty=3, col="grey")
        graphics::axis(side=1, at=Lambda[xseq], labels=formatC(Lambda[xseq], format="e", digits=2), las=2)
        graphics::axis(side=3, at=Lambda[xseq], las=2,
                       labels=formatC(Q[xseq], format="f", big.mark=",", digits=0))
        graphics::mtext(text=xlab, side=1, line=5.2, cex=1.5)
        graphics::mtext(text=zlab, side=3, line=3.5, cex=1.5)
      }
    }

    if (metric=="pi"){
      if (is.null(mar)){
        mar=c(7,5,2,2)
      }
      if (is.null(mfrow)){
        mfrow=c(1,nblocks)
      }
      graphics::par(mar=mar, mfrow=mfrow)

      for (b in block_id){
        # Extracting the stability scores
        if (length(stability$params$pk)==1){
          mat=stability$S_2d
          ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
          mat=mat[ids,,drop=FALSE]
        } else {
          mat=stability$S_2d[,,b]
          ids=which(apply(mat,1,FUN=function(x){any(!is.na(x))}))
          mat=mat[ids,,drop=FALSE]
        }

        # Extracting the best stability score (with optimal lambda) for each pi value
        vect=apply(mat,2,max,na.rm=TRUE)

        # Using input ylab if not as default
        ylab=ifelse(as.character(ylab)==as.character(expression(pi)), yes="Stability Score", no=ylab)

        # Using input xlab if not as default
        xlab=ifelse(as.character(xlab)==as.character(expression(lambda)), yes=expression(pi), no=ylab)

        # Making plot
        cex_points=0.7
        plot(1:length(vect), vect, pch=19, col="navy", cex=cex_points,
             xlab="", ylab=ylab, cex.lab=1.5, xaxt="n")
        xticks=graphics::axTicks(side=1)
        if (min(xticks)==0){
          xticks=xticks+1
        }
        graphics::abline(v=xticks, lty=3, col="grey")
        graphics::abline(h=graphics::axTicks(side=2), lty=3, col="grey")
        graphics::abline(h=max(vect), lty=2, col="red")
        graphics::abline(v=which.max(vect), lty=2, col="red")
        graphics::points(1:length(vect), vect, pch=19, col="navy", cex=cex_points)
        if (lines){
          graphics::lines(1:length(vect), vect, col="navy")
        }

        # Adding x-axis and its labels
        graphics::axis(side=1, at=xticks, labels=formatC(stability$params$pi_list[xticks], digits=2), las=2)
        graphics::mtext(text=xlab, side=1, line=5.2, cex=1.5)
      }
    }
  }

  if (!is.null(filename)){
    grDevices::dev.off()
  }
}


Heatmap=function(mat, colours=c("ivory", "navajowhite", "tomato", "darkred"), resolution=10000, legend=TRUE, legend_length=15, legend_range=NULL){
  # Preparing colours
  colours=grDevices::colorRampPalette(colours)(resolution)
  names(colours)=1:resolution

  # Re-formatting matrix
  mat=mat[,ncol(mat):1]
  vect=as.vector(mat)

  # Defining extreme values
  if (is.null(legend_range)){
    myrange=c(min(vect, na.rm=TRUE), max(vect, na.rm=TRUE))
  } else {
    myrange=legend_range
  }

  # Getting corresponding colours
  mycol=as.character(cut(vect, breaks=seq(myrange[1], myrange[2], length.out=resolution+1), labels=1:resolution, include.lowest=TRUE))
  mycol_mat=matrix(mycol, ncol=ncol(mat))

  # Making heatmap
  plot(NA, xlim=c(0,nrow(mycol_mat)), ylim=c(0,ncol(mycol_mat)),
       xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  for (i in 0:(nrow(mycol_mat)-1)){
    for (j in 0:(ncol(mycol_mat)-1)){
      graphics::polygon(x=c(i,i+1,i+1,i), y=c(j,j,j+1,j+1),
                        col=colours[mycol_mat[i+1,j+1]],
                        border=colours[mycol_mat[i+1,j+1]])
    }
  }

  # Adding colour bar (legend)
  if (legend){
    graphics::par(xpd=TRUE)
    legend_width_factor=1.05
    myrange=round(myrange)
    mylegend_values=seq(myrange[1],myrange[2],length.out=100)
    mylegend_values=unique(round(mylegend_values, digits=-max(nchar(round(myrange)))+2))
    if (is.null(legend_range)){
      mylegend_values=mylegend_values[mylegend_values>=min(myrange)]
    }
    if (length(mylegend_values)>legend_length){
      mylegend_values=unique(round(mylegend_values, digits=-max(nchar(round(myrange)))+1))
    }
    mylegend_ids=as.numeric(as.character(cut(mylegend_values, breaks=seq(myrange[1], myrange[2], length.out=resolution+1),
                                             labels=1:resolution, include.lowest=TRUE)))
    ypos=ncol(mat)
    xpos=nrow(mat)*1.05
    for (l in 1:length(colours)){
      graphics::polygon(x=c(xpos,xpos*legend_width_factor,xpos*legend_width_factor,xpos),
                        y=c(ypos-legend_length+legend_length*l/length(colours),
                            ypos-legend_length+legend_length*l/length(colours),
                            ypos-legend_length+legend_length*(l+1)/length(colours),
                            ypos-legend_length+legend_length*(l+1)/length(colours)),
                        col=colours[l], border=colours[l])
      if (l%in%mylegend_ids){
        graphics::text(x=xpos*legend_width_factor, y=ypos-legend_length+legend_length*(l+0.5)/length(colours),
                       labels=paste0("- ", mylegend_values[which(mylegend_ids==l)]), adj=c(0,0.5))
      }
    }
    graphics::par(xpd=FALSE)
  }
}



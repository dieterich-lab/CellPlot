# EXAMPLE CODE
if (0) {
  x = sort( runif(16, min = 1, max = 3), decreasing = T )
  x.ann = round(runif(16, min=101,max=200))    # THIS IS NEW
  names(x) = paste("GO Term",1:16)
  
  ## Label colors
  xcolor = c(rep("darkslategrey",4), rep("chartreuse4",4), rep("coral4",8))
  
  ## Generate list with random vectors, one for each entry in x
  cells = list()
  xc = round( runif(16, min=21, max=100) )
  for (i in 1:length(xc)) { cells = c(cells, list(runif(xc[i],-5,5))) }
  cells[[9]][1:2] = Inf
  cells[[9]][3] = -Inf
  
  # HERE WE GO
  rs.sym4go( x, x.ann, cells=cells, x.col = xcolor, sort=F, x.mar=c(0,0) )
}


# SYMPLOT FUNCTIONS -- STILL RATHER MESSY

rs.sym4go = function( x, x.annotated, x.up=NULL, x.down=NULL, cells=NULL, x.col=NULL, as.percentage=T, sort=T,
                      x.mar = c(0,0) ) {
  
  # parameter checks
  if ( (is.null(x.down) || is.null(x.up)) && is.null(cells) ) { stop("You must provide either the cells parameter (expression values for members of each term) or numbers of up- and downregulated genes!") }
  if ( !is.null(x.down) ) { if ( (length(x.down) != length(x.up) ) || (length(x.down) != length(x) ) ) { stop("x, x.up, and x.down must contain equal numbers of elements.")} }
  if ( length(x) != length(x.annotated) ) {  stop("x and x.annotated must contain an equal number of elements.") }
  
  
  if (!is.null(cells) ) {
    x.up = sapply(cells, function(x) sum(x > 0) )
    x.down = sapply(cells, function(x) sum(x < 0) )
  }
  
  outframe = matrix(0, nrow=length(x), ncol=5, dimnames = list(names(x),NULL) )
    for (i in 1:length(x)) {
      outframe[i,] = c(x.up[i],0,x.annotated[i],0,x.down[i])
    }
    if (as.percentage) {
      #outframe = t( apply(outframe,1,function(x) { x[c(1,5)] = x[c(1,5)] / x[3] * 100 ; return(x) } ) )
    }
  outframe = as.data.frame(outframe)
  outframe$val = x
  if (!is.null(x.col)) { outframe$col = x.col }
 
  if (sort) {
    outframe = outframe[order(outframe[,3]),]
  } else {
    outframe = outframe[nrow(outframe):1,]
  }
  
  rs.symplot( outframe[,1:5], label.col = outframe[,7], ticksize = 10, group.labels = c("Downregulated","Annotated","Upregulated"), inclusion = T, gridlines = T,
              x.mar=x.mar, mid.values=outframe[,6], as.percentage=T, mid.bounds=c(0,ceiling(max(outframe[,6]))), key.lab="log2(GO term enrichment)", bar.scale=1.4)
}


rs.symplot = function( symframe, x.mar=c(0.2,0.05), y.mar = c(0.1,0), bar.lwd=2, bar.scale=1, shading.density=20, space = 0.1, midgap=0.1, mid.values=NULL, 
                       mid.bounds=NULL, mid.col=c("white","darkred"), label.col=NULL, key.lab="color key", label.align.left=T, gridlines=T,
                       cols = c("deepskyblue2","grey","coral"), group.labels=c("Downregulated","Annotated","Upregulated"), 
                       group.cex=0.8, axis.cex=0.8, mid.cex=0.8, label.cex=0.8, ticksize=NULL, xlim=NULL, inclusion=F, as.percentage=F, ...) {
  
  if(!is.null(mid.bounds)) {
    midbound=mid.bounds
  } else {
    midbound = c( min(mid.values), max(mid.values) )
  }
  
  if(!is.null(mid.values)) {
    midcolmap = seq( midbound[1], midbound[2], length.out=101 )
    names(midcolmap) = colorRampPalette( c(mid.col[1], mid.col[2]) )(101)
    mid.i = sapply(mid.values, function(x) which.min(x > midcolmap))
    mid.c = names(midcolmap)[mid.i]
  } else {
    mid.c = rep( mid.col[1], nrow(symframe) )
  }
  
  if (as.percentage) {
    symframe = t(apply(symframe,1,function(x) {
      x[1:2] = x[1:2] / x[3] *100
      x[4:5] = x[4:5] / x[3] *100
      return(x)
    }))
    xlim=max(symframe[,c(1,2,4,5)])
    xlim=xlim+xlim%%ticksize
  }
  
  par(xpd=NA)
  if(is.null(xlim)) { 
    bound = max( apply( symframe,1,function(x) max( sum(x[1:2]), sum(x[4:5])) ) )
    bound = bound + bound %% ticksize 
  } else {
    bound = xlim
  }
  boundmid = max(symframe[,3])
  midmean = boundmid/2
  
  ybound = c(1,0) + c(-1,1)*y.mar
  yscale = par("pin")[1]/par("pin")[2] * diff(par("usr")[1:2])/diff(par("usr")[3:4])
  ybound[2] = ybound[1] - ( (ybound[1] - ybound[2]) * yscale * 0.4 * bar.scale )
  
  # if (!is.null(scaleTo)) { ybound[2] = ybound[1] - ( (ybound[1] - ybound[2]) / scaleTo ) * nrow(symframe) }
  xbound = c(0,1) + c(1,-1)*x.mar
  x.left = c( xbound[1], (xbound[2]-xbound[1])*(0.5-midgap)+xbound[1] )
  x.right = c( (xbound[2]-xbound[1])*(0.5+midgap)+xbound[1], xbound[2] )
  ysteps = seq( ybound[2], ybound[1], length.out=( nrow(symframe)+1 ) )
  
  ygap = abs(ysteps[1]-ysteps[2])
  yspace = space * ygap
  
  
  if (!is.null(ticksize)) {
    ticklabels = round( seq(0,bound,ticksize), digits = 0 )
    left.axis.at = (sort(ticklabels,decreasing = T)/bound) * (x.left[2]-x.left[1]) + x.left[1] + (bound-max(ticklabels))/bound * (x.left[2]-x.left[1])
    right.axis.at = (ticklabels/bound) * (x.right[2]-x.right[1]) + x.right[1]
  } else {
    ticklabels = round( seq(0,bound,bound/ticks), digits = ifelse(scale,1,0) )
    left.axis.at = seq(x.left[2], x.left[1], by= -( (x.left[2] - x.left[1]) / ticks)  )
    right.axis.at = seq(x.right[1], x.right[2], by= ( (x.right[2] - x.right[1]) / ticks)  )
  }
  if (as.percentage) { ticklabels = paste0(ticklabels, "%")}
  
  
  # BARS
  plot.new()
  
  # GRID
  if (gridlines) {
    segments(c(left.axis.at,right.axis.at), ybound[2]-yspace, c(left.axis.at,right.axis.at), ybound[1]+yspace, col="grey", lwd = bar.lwd )
    segments(left.axis.at[length(left.axis.at)],ybound[2]-yspace,left.axis.at[1],ybound[2]-yspace, col="grey", lwd = bar.lwd)
    segments(right.axis.at[length(right.axis.at)],ybound[2]-yspace,right.axis.at[1],ybound[2]-yspace, col="grey", lwd = bar.lwd)
  }
  
  
  
  for (i in 1:nrow(symframe)) {
    
    # the one in the middle
    midspace = 0.1 * (x.right[1] - x.left[2])
    midrange = (x.right[1] - x.left[2]) - 2 * midspace
    #midspace = (midspace + (boundmid - symframe[i,3])) * midrange/ boundmid + midrange*0.1
    midspace = ((boundmid - symframe[i,3]) / (2*boundmid) ) * midrange + midspace
    
    
    rect( x.left[2], ysteps[i+1] - yspace, x.right[1], ysteps[i] + yspace, border = NA, col=cols[2], density = shading.density, lwd = 2 )
    # rect( x.left[2]+midspace, ysteps[i+1] - yspace, x.right[1]-midspace, ysteps[i] + yspace, col="darkgray", lwd = bar.lwd, density = shading.density)
    rect( x.left[2]+midspace, ysteps[i+1] - yspace, x.right[1]-midspace, ysteps[i] + yspace, lwd = bar.lwd, col=mid.c[i])
    
    # pillars of the community
    #segments(x.left[2], ybound[2], x.left[2], ybound[1])
    #segments(x.right[1], ybound[2], x.right[1], ybound[1])
    
    # left all
    rect( (bound-sum(symframe[i,1:2]))/bound * (x.left[2]-x.left[1]) + x.left[1], 
          ysteps[i+1] - yspace,
          x.left[2],
          ysteps[i] + yspace,
          col=cols[1], lwd = bar.lwd)
    # left subset
    rect( (bound-symframe[i,2])/bound * (x.left[2]-x.left[1]) + x.left[1], 
          ysteps[i+1] - yspace,
          x.left[2],
          ysteps[i] + yspace,
          col = "black", lwd = bar.lwd, density=shading.density)
    
    # right all
    rect( x.right[1],
          ysteps[i+1] - yspace,
          sum(symframe[i,4:5])/bound*(x.right[2]-x.right[1]) + x.right[1], 
          ysteps[i] + yspace,
          col=cols[3], lwd = bar.lwd)
    # right subset
    rect( x.right[1],
          ysteps[i+1] - yspace,
          (sum(symframe[i,4:5])-symframe[i,5]) /bound*(x.right[2]-x.right[1]) + x.right[1], 
          ysteps[i] + yspace,
          col = "black", lwd = bar.lwd, density=shading.density)
    
    # TEXT   
    # text middle
    if (symframe[i,3] > midmean ) {
      text( (xbound[2]-xbound[1])/2 + xbound[1], ysteps[i]+0.5*ygap, symframe[i,3], 
            col = "black", font=2, cex = mid.cex )
    } else { 
      text( x.right[1]-midspace+0.09*midrange, ysteps[i]+0.5*ygap, symframe[i,3], 
            col = "black", font=2, cex = mid.cex ) 
    }
    
    
    # labels
    text( x.left[1]+(x.left[2]-x.left[1])*0.05, ysteps[i+1]-0.5*ygap, pos=2, rownames(symframe)[i], cex = label.cex, col=ifelse(is.null(label.col),"black",label.col[i] ))
  }
  
  # I AM LEGEND
  ylegend = ybound[1] + 0.04*yscale
  
  rect(x.left[1]+(x.left[2]-x.left[1])*0.1, ylegend+ygap*1.8, x.left[2], ylegend+ygap*0.7, col=cols[1], lwd = bar.lwd)
  rect(x.right[1], ylegend+ygap*1.8, x.right[2]-(x.right[2]-x.right[1])*0.1, ylegend+ygap*0.7, col=cols[3], lwd = bar.lwd)
  
  if (!inclusion) {
    #rect(x.left[2]-midrange*0.5, ylegend+ygap*2, x.right[1]+midrange*0.5, ylegend+ygap*0.5, lwd = bar.lwd)
    rect(x.left[2]-(x.left[2]-x.left[1])*0.2, ylegend+ygap*2, x.right[1]+(x.left[2]-x.left[1])*0.2, ylegend+ygap*0.5, lwd = bar.lwd)
    
    rect(x.left[2]-(x.left[2]-x.left[1])*0.2, ylegend+ygap*1.8, x.left[2], ylegend+ygap*0.7, lwd = bar.lwd, col="black", density=shading.density)
    rect(x.right[1], ylegend+ygap*1.8, x.right[1]+(x.left[2]-x.left[1])*0.2, ylegend+ygap*0.7, lwd = bar.lwd, col="black", density=shading.density)
    
  } else {
    rect(x.left[1]+(x.left[2]-x.left[1])*0.05, ylegend+ygap*2, x.right[2]-(x.right[2]-x.right[1])*0.05, ylegend+ygap*0.5, lwd = bar.lwd)
  }
  
  text( x.left[1]+(x.left[2]-x.left[1])*0.1, ylegend + ygap * 1.25, group.labels[1], pos = 4, font=2, cex = group.cex )
  text( x.right[2]-(x.right[2]-x.right[1])*0.1, ylegend + ygap * 1.25, group.labels[3], pos = 2, font=2, cex = group.cex )
  text( (xbound[2]-xbound[1])/2 + xbound[1], ylegend + ygap * 1.25, group.labels[2], font=2, cex = group.cex )
  
  # AXIS
  axis(3, pos=ybound[1]+yspace, at = left.axis.at, labels = ticklabels, cex.axis=axis.cex, lwd = bar.lwd, padj = 1 )
  axis(3, pos=ybound[1]+yspace, at = right.axis.at, labels = ticklabels, cex.axis=axis.cex, lwd = bar.lwd, padj = 1 )  
  
  # TITLE
  title(...)
  
  # COLOR LEGEND IF MID VALUES ARE PROVIDED
  # color legend
  key.n=11
  if (!is.null(mid.values)) {
    #midbound = range(mid.values)
    lc = c( left.axis.at[length(left.axis.at)],ybound[2]-2*yspace-ygap,right.axis.at[length(right.axis.at)],ybound[2]-4*yspace )
    absmax = max(abs(midbound))
    lc.min <- min(midbound)
    lc.max <- max(midbound)
    midcolmap.center <- midcolmap[51]
    midcolmap <- midcolmap[lc.min <= midcolmap & midcolmap <= lc.max]
    lc.xsteps = seq( lc[1], lc[3], length.out=key.n+1 )
    lc.xgap = lc.xsteps[1] - lc.xsteps[2]
    
    lc.density <- rep(0, key.n)
    lc.range = seq( min(lc.min,midbound[1]), lc.max, length.out=key.n )
    lc.col <- names(midcolmap)[seq(1,length(midcolmap),length.out=key.n)]
    
    # normalize center zero color and text
    if (key.n %% 2 > 0 && lc.min < 0 && lc.max > 0) {
      i.center <- mean(c(1,key.n))
      lc.col[i.center] <- names(midcolmap.center)
      lc.range[i.center] <- 0
    }
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col=lc.col, lwd=bar.lwd )
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col="black", lwd=bar.lwd, density = lc.density, border=NA )
    text( (lc.xsteps[-(key.n+1)]+lc.xsteps[-1])/2, lc[2], pos=1, labels=round(lc.range,1), cex=mid.cex, font=2 )
    text( (xbound[1]+xbound[2])/2, lc[2]-strheight("0",cex=mid.cex)*1.5 , labels=key.lab, pos=1, cex=mid.cex )
  }
}
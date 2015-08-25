#' @title Sym Plot
#'
#' @description 
#' Plots a split barchart, showing the proportions of two mutually exclusive
#' sets in relation to a set containing them both. E.g., Gene Ontology terms,
#' showing the proportions of differentially down-regulated and up-regulated
#' annotated genes from a perturbation experiment.
#' The color of the central column elements maps to the value provided in x
#' (e.g. GO term enrichment). Associated genes may be provided as a list of
#' vectors of expression values, same as for cell.plot(), or as separate
#' vectors x.up and x.down, providing the numbers of up- and down-regulated
#' genes in the same order as x.
#'
#' @param x Strictly positive numeric vector.
#'
#' @param x.annotated Cardinality of each functional set (e.g. total number of genes annotated
#' in a given GO term). Must be the same length and order as x.
#'
#' @param x.up Vector containing numbers of up-regulated genes in each functional set. Must be the
#' same length and order as x. Overridden by the cells parameter, if provided.
#' 
#' @param x.down Vector containing numbers of down-regulated genes in each functional set. See x.up.
#'
#' @param cells List of vectors (e.g. gene logfold-changes). Must be the same length and order as x.
#' 
#' @param sort Sort categories by their cardinality. Defaults to FALSE.
#' 
#' @param main Title.
#' 
#' @param elem.bounds Vector of length 2 specifying a filter on minimum and maximum number of elements
#' in both categories combined. Defaults to NULL, in which case no filter is applied.
#'
#' @param x.mar Left and right margins as fractions of plot width. Defaults to c(0.2,0.05).
#'
#' @param y.mar Lower and upper margins as fractions of plot height. Defaults to c(0.1,0).
#'
#' @param bar.lwd Line width of bar elements. Defaults to 2.
#'
#' @param bar.scale If not NULL (the default), set bar height to a fixed value relative to the
#' numeric factor provided. Use to ensure consistency across plots with varying numbers of elements.
#' 
#' @param space Free space above and below bars as a fraction of bar height. Defaults to 0.1.
#' 
#' @param mid.gap Free space left between the extremes of the central column and the bar elements, as a
#' fraction of central column width. Defaults to 0.1.
#' 
#' @param mid.bounds Lower and upper bound on the color scale mapping of the central column 
#' (i.e. mapping of x). Defaults to c(0, max(x)).
#' 
#' @param mid.col Color scale for mapping of x to the central column elements. Defaults to 
#' c("white","darkred").
#' 
#' @param key.lab Label for color key.
#'
#' @param key.n Number of legend boxes for the color key. Defaults to 11.
#'
#' @param cols Colors for left and right bars. Defaults to c("deepskyblue2","coral").
#' 
#' @param group.labels Labels for left bars, central column, and right bars. Defaults to
#' c("Downregulated","Annotated","Upregulated").
#'
#' @param group.cex Scaling factor for group labels. Defaults to 0.8.
#' 
#' @param axis.cex Scaling factor for axis labeling. Defaults to 0.8.
#' 
#' @param mid.cex Scaling factor for central column labeling. Defaults to 0.8.
#' 
#' @param lab.cex Scaling factor for functional group labels. Defaults to 1.
#' 
#' @param ticksize Spacing between x-axis ticks. Defaults to 10.
#'
#' @author 
#' Robert Sehlke [aut]\cr
#' Sven E. Templer [ctb]
#'
#' @examples
#' \dontrun{
#' ## Generate random positive vector and name it
#' x = sort( runif(16, min = 1, max = 3), decreasing = T )
#' names(x) = paste("GO Term",1:16)
#'
#' ## Label colors
#' xcolor = c(rep("darkslategrey",4), rep("chartreuse4",4), rep("coral4",8))
#'
#' ## Generate list with random vectors, one for each entry in x
#' cells = list()
#' xc = round( runif(16, min=21, max=100) )
#' for (i in 1:length(xc)) { cells = c(cells, list(runif(xc[i],-5,5))) }
#' cells[[9]][1:2] = Inf
#' cells[[9]][3] = -Inf
#' cells.annotated <- sapply(cells, length)
#'
#' ## Plot with spacers
#' sym.plot( x, cells, cells.annotated, 
#'           xcolor, spacers = c(4,8), xlab.ticks = 5, cell.limit = 80,
#'           main="Cell Plot Demo", xlab="log(enrichment)" )
#'   
#' ## golub.deg data example:
#' data(golub.deg)
#' sym.plot( x = golub.deg$go$go.loge, cells = golub.deg$go$deg.logfc,
#'           x.annotated = golub.deg$go$go.table$Annotated, 
#'           elem.bounds = c(5,100), bar.scale=1, x.mar=c(0.3,0),
#'           main = "Golub et al. DEG and GO Enrichment" )
#' }
#' 
#'

#' @export
sym.plot = function( x, cells=NULL, x.annotated, x.up=NULL, x.down=NULL, x.col=NULL, sort=F, main="", elem.bounds=NULL, 
                     x.mar=c(0.2,0.05), y.mar = c(0.1,0), bar.lwd=2, bar.scale=NULL, space = 0.1, mid.gap=0.1,
                     mid.bounds=NULL, mid.col=c("white","darkred"), key.lab="GO Term Enrichment", key.n = 11, cols = c("deepskyblue2","coral"), 
                     group.labels=c("Downregulated","Annotated","Upregulated"), group.cex=0.8, axis.cex=0.8, mid.cex=0.8, 
                     lab.cex=1, ticksize=10, xlim=NULL, ...) {
  shading.density=20
  
  # parameter checks
  if ( (is.null(x.down) || is.null(x.up)) && is.null(cells) ) { stop("You must provide either the cells parameter (expression values for members of each term) or numbers of up- and downregulated genes!") }
  if ( !is.null(x.down) ) { if ( (length(x.down) != length(x.up) ) || (length(x.down) != length(x) ) ) { stop("x, x.up, and x.down must contain equal numbers of elements.")} }
  if ( length(x) != length(x.annotated) ) {  stop("x and x.annotated must contain an equal number of elements.") }
  
  
  if (!is.null(cells) ) {
    x.up = sapply(cells, function(x) sum(x > 0) )
    x.down = sapply(cells, function(x) sum(x < 0) )
  }
  
  # formatting input
  if(!is.null(elem.bounds)) {
    ec = x.up + x.down
    excl = which( (ec < elem.bounds[1]) | (ec > elem.bounds[2]) )
    if (length(excl) == length(x) ) { stop("No elements in the specified range!") }
    if (length(excl) > 0) {
      x = x[-excl]
      x.annotated = x.annotated[-excl]
      x.up=x.up[-excl]
      x.down=x.down[-excl]
      if (!is.null(x.col)) { x.col = x.col[-excl] }
    }
  }
  
  outframe = matrix(0, nrow=length(x), ncol=5, dimnames = list(names(x),NULL) )
  for (i in 1:length(x)) {
    outframe[i,] = c(x.down[i],0,x.annotated[i],0,x.up[i])
  }
  outframe = as.data.frame(outframe)
  outframe$val = x
  if (!is.null(x.col)) { outframe$col = x.col } else { outframe$col = rep("black",nrow(outframe)) }
  
  if (sort) {
    outframe = outframe[order(outframe[,3]),]
  } else {
    outframe = outframe[nrow(outframe):1,]
  }
  
  if(is.null(mid.bounds)) { mid.bounds = c(0,ceiling(max(outframe[,6]))) }
  symframe = outframe[,1:5]
  label.col = outframe[,7]
  inclusion = T
  as.percentage=T
  gridlines = T
  mid.values=outframe[,6]
  yscale = (diff(par("usr")[3:4])/par("pin")[2])
  
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
    xl=max(symframe[,c(1,2,4,5)])
    while ( xl/ticksize < 1 ) {
      ticksize=ticksize/10
    }
    xlim=min( xl+ticksize%%xl, 100)
  }
  
  par(xpd=NA)
  if(is.null(xlim)) { 
    bound = max( apply( symframe,1,function(x) max( sum(x[1:2]), sum(x[4:5])) ) )
    bound = bound + bound %% ticksize 
  } else {
    bound = xlim
  }
  boundmid = max(symframe[,3])
  midmean = boundmid/3
  
  ybound = c(1,0) + c(-1,1)*y.mar
  xbound = c(0,1) + c(1,-1)*x.mar
  x.left = c( xbound[1], (xbound[2]-xbound[1])*(0.5-mid.gap)+xbound[1] )
  x.right = c( (xbound[2]-xbound[1])*(0.5+mid.gap)+xbound[1], xbound[2] )
  ysteps = seq( ybound[2], ybound[1], length.out=( nrow(symframe)+1 ) )
  
  if ( !is.null(bar.scale) ) {
    ybound[2] = ybound[1]-( 0.3 * yscale * bar.scale * nrow(symframe) )
    ysteps = ybound[2] + c(0, cumsum( rep(0.3 * yscale * bar.scale, nrow(symframe) ) ) )
    
    if ( ybound[2] < par("usr")[3] ) {
      warning("Plotting area too small! Decrease bar.scale or increase vertical space.")
    }
  }
  
  ygap = abs(ysteps[1]-ysteps[2])
  yspace = space * ygap
  
    ticklabels = round( seq(0,bound,ticksize), digits = ifelse(ticksize>=1, 0, 1) )
    left.axis.at = (sort(ticklabels,decreasing = T)/bound) * (x.left[2]-x.left[1]) + x.left[1] + (bound-max(ticklabels))/bound * (x.left[2]-x.left[1])
    right.axis.at = (ticklabels/bound) * (x.right[2]-x.right[1]) + x.right[1]

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
    midspace = ((boundmid - symframe[i,3]) / (2*boundmid) ) * midrange + midspace
    
    rect( x.left[2], ysteps[i+1] - yspace, x.right[1], ysteps[i] + yspace, border = NA, col="grey", density = shading.density, lwd = 2 )
    rect( x.left[2]+midspace, ysteps[i+1] - yspace, x.right[1]-midspace, ysteps[i] + yspace, lwd = bar.lwd, col=mid.c[i], border = mid.c[i])
    
    # left all
    rect( (bound-sum(symframe[i,1:2]))/bound * (x.left[2]-x.left[1]) + x.left[1], 
          ysteps[i+1] - yspace,
          x.left[2],
          ysteps[i] + yspace,
          col=cols[1], border=cols[1], lwd = bar.lwd)
    # left subset
#     rect( (bound-symframe[i,2])/bound * (x.left[2]-x.left[1]) + x.left[1], 
#           ysteps[i+1] - yspace,
#           x.left[2],
#           ysteps[i] + yspace,
#           col = "black", lwd = bar.lwd, density=shading.density)
    
    # right all
    rect( x.right[1],
          ysteps[i+1] - yspace,
          sum(symframe[i,4:5])/bound*(x.right[2]-x.right[1]) + x.right[1], 
          ysteps[i] + yspace,
          col=cols[2], border=cols[2], lwd = bar.lwd)
    # right subset
#     rect( x.right[1],
#           ysteps[i+1] - yspace,
#           (sum(symframe[i,4:5])-symframe[i,5]) /bound*(x.right[2]-x.right[1]) + x.right[1], 
#           ysteps[i] + yspace,
#           col = "black", lwd = bar.lwd, density=shading.density)
    
    # TEXT   
    # text middle
    if (symframe[i,3] > midmean ) {
      text( (xbound[2]-xbound[1])/2 + xbound[1], ysteps[i]+0.5*ygap, symframe[i,3], 
            col = "black", font=2, cex = mid.cex )
    } else { 
      text( x.right[1]-0.3*midrange, ysteps[i]+0.5*ygap, symframe[i,3], 
            col = "black", font=2, cex = mid.cex ) 
    }
    
    
    # labels
    text( x.left[1]+(x.left[2]-x.left[1])*0.05, ysteps[i+1]-0.5*ygap, pos=2, rownames(symframe)[i], cex = lab.cex, col=ifelse(is.null(label.col),"black",label.col[i] ))
  }
  
  # I AM LEGEND
  ylegend = ybound[1] + 0.4*yscale
  
  rect(x.left[1]+(x.left[2]-x.left[1])*0.1, ylegend+ygap*1.8, x.left[2], ylegend+ygap*0.7, col=cols[1], border=cols[1], lwd = bar.lwd)
  rect(x.right[1], ylegend+ygap*1.8, x.right[2]-(x.right[2]-x.right[1])*0.1, ylegend+ygap*0.7, col=cols[2], border=cols[2], lwd = bar.lwd)
  
  if (!inclusion) {
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
  title(main=main)
  
  # COLOR LEGEND IF MID VALUES ARE PROVIDED
  # color legend
  if (!is.null(mid.values)) {
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
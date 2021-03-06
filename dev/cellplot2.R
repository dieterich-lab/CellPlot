#' @name cellplot2
#' @title Cell Plot (2)
#'
#' @description 
#' Creates a barplot with bars containing a heatmap of variable amount of 
#' cells. Can be used for example to display the log enrichment of gene 
#' ontology (GO) terms in a subset of differentially regulated genes.
#' The heatmap then can present the log fold change of all significantly
#' regulated genes within the respective GO term.
#' 
#' @param x Bar size vector (e.g. log of GO enrichment).
#' @param classes Names of the classes (e.g. GO terms).
#' @param cells List of numeric vectors of same size as x.

##### subset.cell.n.range
#' @param elem.bounds Vector of length 2 specifying a filter on minimum and maximum number of elements
#' in a category. Defaults to NULL, in which case no filter is applied. 

##### xlim2
#' @param x.bound Numeric, must be positive. Specifies upper bound of x-axis scale. 
#' If NULL (the default), this value is chosen automatically.

#' @param bar.scale Numeric. Set bar height to a fixed value multiplied by this
#' parameter. If the plotting area is too small clipping may occur.

#' @param cell.col Character vector of length 3, specifying colors for the low 
#' end, zero, and high end of cell values. Defaults to c("blue","white","red").
#' @param inf.shading Numeric. Density of shading lines for infinite values. 
#' Defaults to 30/cell.lwd.

#' @param space Scaling factor for spacing between bars.
#' @param x.mar Numeric vector defining the inside margin to the left and right
#' of the main plot. Defaults to c(0.2,0.1).
#' @param y.mar Numeric vector defining the inside margin to the top and bottom
#' of the main plot. Defaults to c(0.08,0).


#' @param lab.cex Scaling factor for label text size.
#' @param xdes.cex Scaling factor for key legend text size.
#' @param xlab.cex Scaling factor for x-label text size.
#' @param xlab.ticks Number of ticks for the x-axis..

#' @param cell.lwd Size of the border of individual cells. Defaults to 1.
#' @param cell.outer Size of outer border around bars. Defaults to 2.
#' @param cell.sort Should cell values be sorted? Defaults to TRUE.
#' @param cell.limit Number of cells above which separators between cells are
#' omitted. Defaults to 50.
#' @param cell.bounds Numeric vector of length 2, specifying lower and upper 
#' bound for cell value visualization. If set to NULL (the default), these
#' values are chosen automatically.

#' @param xlab Label for x-axis.
#' @param key Logical value. Show the key/legend for the cell colors.
#' @param key.lab Label for color key.
#' @param key.n Number of legend boxes for the color key. Minimum is 3 (5 when
#' infinite values are present).
#' @param spacers Numeric vector. Inserts empty space after the specified
#' positions to visually group bars.
#' @param sym Logical, if \code{TRUE} visualize cell values on a symmetrical scale.
#' @param gridlines Logical value, if to show vertical grid lines.
#' @param \dots Arguments passed through to the title() function. E.g.
#' main = "Plot title".
#'
#' @author 
#' Robert Sehlke [aut]\cr
#' Sven E. Templer [ctb]


#' @rdname cellplot2
#' @export
cellplot2 <- function (height, cells, ...) 
  UseMethod("cellplot2")

#' @rdname cellplot2
#' @export
cellplot2.data.frame <- function (height, cells, ...)
  stop("not yet developed")
  
#' @rdname cellplot2
#' @export
cellplot2.integer <- function (height, cells, ...)
{
  x <- as.numeric(x)
  NextMethod("cellplot2")
}


#' @rdname cellplot2
#' @export
cellplot2.numeric <- function(
  height, cells, class.names = NULL,
  cells.col = c("deepskyblue2", "white", "coral"),
  cells.col.centered = TRUE,
  cells.sort = TRUE, cells.sort.decreasing = FALSE,
  horiz = TRUE, width = .8,
  legend.n.keys = 5,
  legend.pos = "right", # see ?xy.coords, e.g. matrix(0:1,,2)
  legend.n.sign = 5,
  barplot.mar = c(.1,.1,.1,.1),
  axis.n.ticks = 5,
  
  x, cells, x.col=NULL, cell.col=,
  inf.shading = 30/cell.lwd,  space=0.1, x.mar=c(0.2, 0), y.mar = c(0.1, 0), x.bound=NULL, lab.cex = 1, xdes.cex=1, xlab.cex=1, xlab.ticks=5,
  sym=FALSE, cell.lwd=2, cell.outer=2, cell.sort=T, cell.limit=30, cell.bounds=NULL, elem.bounds=NULL, xlab="GO Term Enrichment",
  key=T, key.lab="Differential Expression", key.n=11, spacers=NULL, bar.scale=1, gridlines=T, ... )
{
  
  ### sampling
  height <- sort(sample(4:11,5), decreasing = F)
  cells <- lapply(height, function(i) rnorm(sample(3:10,1)))
  horiz <- TRUE
  width <- .8
  cells.col.centered <- TRUE
  cells.sort.decreasing = FALSE
  cells.sort = TRUE
  cells.col <- c("deepskyblue2", "white", "coral")
  legend.n.keys = 5
  legend.pos = "right"
  legend.n.sign = 5
  barplot.mar = c(0,.2,0,.1)
  axis.n.ticks = 5
  ### end sampling
  
  # attributes
  n <- length(height)
  hmax <- max(height)
  hr <- range(c(0, height))
  pos <- seq(n)      # bar positions
  width <- width[1]/2
  
  # convert
  cells <- lapply(cells, as.numeric)
  if (is.null(class.names))
    class.names <- names(height)
  if (!is.null(class.names))
    class.names <- as.character(class.names)
  if (is.null(class.names))
    class.names <- as.character(seq.int(n))
  if (cells.sort)
    cells <- lapply(cells, sort, decreasing = cells.sort.decreasing)
  
  # attributes
  cn <- sapply(cells, length)
  cr <- range(unlist(cells), na.rm = TRUE)
  
  # checks
  if (any(is.na(height)))
    stop("missing values in height")
  if (!is.list(cells))
    stop("cells is not a list")
  if (length(cells) != n)
    stop("length of cells not equal to x")
  if (length(class.names) != n)
    stop("length of class.names not equal to x")
  if (any(cn < 1))
    stop("cells must at least be of length 1 each")
  if (!cp_within(width, c(0,1)))
    stop("width must be within 0 and 1")
  
  #pos <- rev(as.vector(barplot(h, horiz = horiz)))
  #dev.off()
  plot.new()
  usr <- par('usr')
  usr <- usr + barplot.mar * c(1,-1,1,-1)
#   i <- 5
#   hi <- height[i]
#   pi <- pos[i]
#   ci <- cells[[i]]
#   cp_cells(hi = hi, pi = pi, ci = ci, horiz = horiz, w = width, p = pos, u = usr, h = height)
  null <- Map(cp_cells, 
      hi = height, pi = pos, ci = cells, 
      w = list(width), p = list(pos), u = list(usr), h = list(height), 
      horiz = horiz, cols = list(cells.col), cellrange = list(cr), centered = cells.col.centered)
  
  leg <- if (cells.col.centered) max(abs(cr)) * c(-1, 1) else cr
  leg <- seq(leg[1], leg[2], length.out = legend.n.keys)
  leg <- setNames(leg, cp_colorpal(leg, cols = cells.col, centered = cells.col.centered))
  legend(legend.pos, legend = signif(leg, legend.n.sign), fill = names(leg))
  
  # bottom
  ax.where <- if (horiz) 1 else 3
  ax.usr <- if (horiz) usr[1:2] else usr[3:4]
  ax.val <- seq(hr[1], hr[2], length.out = axis.n.ticks)
  ax.range <- scaler(ax.val, ax.usr)
  axis(ax.where, ax.range, ax.val)
  
  
  
  
  
  
  
  return(NULL)
  
  ###############
  
  
  
  # parameter checks
  if(!is.null(x.bound)){ if(!(is.numeric(x.bound) && (x.bound > 0)) ) {
    stop("x.bound must be a positive numeric value")
  }}
  
  if(!is.null(elem.bounds)) {
    ec = sapply(cells,length)
    excl = which( (ec < elem.bounds[1]) | (ec > elem.bounds[2]) )
    if (length(excl) == length(x) ) { stop("No elements in the specified range!") }
    if (length(excl) > 0) {
      x = x[-excl]
      if (!is.null(cells)) { cells = cells[-excl] } 
      else { x.up=x.up[-excl]; x.down=x.down[-excl] }
      if (!is.null(x.col)) { x.col = x.col[-excl] }
    }
  }
  
  # yscale = par("pin")[1]/par("pin")[2] * diff(par("usr")[1:2])/diff(par("usr")[3:4])
  yscale = (diff(par("usr")[3:4])/par("pin")[2])
  
  par(xpd=NA)
  cell.col.inf = cell.col[c(1,3)]
  #if (is.null(xlab.ticks)) { xlab.ticks = round( max(x) / 10, digits = 1) }
  #ticksize = xlab.ticks
  ybound = c(1,0) + c(-1,1)*y.mar
  
  # scale  
  ysteps = ybound[1] - cumsum( rep(0.3 * yscale * bar.scale, length(x)+1+ifelse(is.null(spacers), 0, length(spacers)) ) )
  ybound[2] = min(ysteps, ybound[2])
  if ( ybound[2] < par("usr")[3] ) {
    warning("Plotting area too small! Decrease bar.scale.fixed or increase vertical space.")
  }
  
  xbound = c(0,1) + c(1,-1)*x.mar
  if (is.null(spacers)) {
    ysteps = seq( ybound[1], ybound[2], length.out=( length(x)+1 ) )
  } else {
    spacers = spacers + 1:length(spacers) + 1
    ysteps = seq( ybound[1], ybound[2], length.out=( length(x)+1+length(spacers) ) )
    ysteps = ysteps[-spacers]
  }
  
  ygap = abs(ysteps[1]-ysteps[2])
  yspace = space * ygap
  
  celldata = unlist(cells)
  cellinf <- is.infinite(celldata)
  cellmis <- is.na(celldata)
  cellbound = range(celldata[!cellinf & !cellmis])
  if (!is.null(cell.bounds) & is.numeric(cell.bounds) ) { cellbound = cell.bounds }
  if (sym) { cellbound = rep( max(abs(cellbound)),2 ) * c(-1,1) }
  # cellcolmap = seq( cellbound[1], cellbound[2], length.out=101 )  
  cellcolmap = c( seq( cellbound[1], 0, length.out=50 ), 0, seq( 0, cellbound[2], length.out=50 ) )
  names(cellcolmap) = c( colorRampPalette( c(cell.col[1], cell.col[2]) )(50),
                         cell.col[2],
                         colorRampPalette( c(cell.col[2], cell.col[3]) )(50) )
  
  if (any(cellinf)) {
    if (key.n < 5) stop("key.n must be >4 when infinite values are present")
  } else {
    if (key.n < 3) stop("key.n must be >2")
  }
  
  
  labvec = rep("",length(x))
  if (!is.null(names(x))) { labvec=names(x) }
  colvec = x.col
  if( is.null(colvec) ) { colvec = rep("black",length(x)) }
  
  # do the actual plotting
  plot.new()
  
  axis.at <- seq(xbound[1], xbound[2], length.out = xlab.ticks)
  axis.lab <- round(seq(min(0,min(x)), max(x), length.out = xlab.ticks),1)
  if (!is.null(x.bound) && is.numeric(x.bound) && (x.bound > 0)) { axis.lab <- round(seq(min(0,min(x)), x.bound, length.out = xlab.ticks),1) }
  
  # GRID
  if (gridlines) {
    segments(axis.at, ybound[2]-yspace, axis.at, ybound[1]+0.015, col="grey", lwd = cell.outer )
    #segments(left.axis.at[length(left.axis.at)],ybound[2]-yspace,left.axis.at[1],ybound[2]-yspace, col="grey", lwd = bar.lwd)
    segments(axis.at[length(axis.at)],ybound[2]-yspace,axis.at[1],ybound[2]-yspace, col="grey", lwd = cell.outer)
  }
  axis(3, pos=ybound[1]+0.015, at = axis.at, labels = axis.lab, cex.axis=xlab.cex, padj=0.5, lwd=cell.outer )
  
  for (i in 1:length(x)) {
    bar.n <- length(cells[[i]])
    bar.nreal <- sum(!is.na(cells[[i]]))
    xsteps = seq(xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1], length.out=(bar.n+1))
    if (!is.null(x.bound) && is.numeric(x.bound) && (x.bound > 0)) { 
      xsteps = seq(xbound[1], (x[i]/x.bound)*(xbound[2]-xbound[1])+xbound[1], length.out=(bar.n+1)) }
    
    if (is.null(cells)) { xsteps = c(xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1]) }
    # row labels
    text( xbound[1], ysteps[i+1]+ygap*0.5, labels=labvec[i], pos=2, cex=lab.cex, col = colvec[i] )
    if(cell.sort) { cells[[i]] = sort(cells[[i]]) }
    bar.order <- order(cells[[i]])
    # number of cells labels
    text( xsteps[length(xsteps)], ysteps[i+1]+ygap*0.5, pos=4, cex=lab.cex, labels=bar.nreal)
    # map colors to cell values
    bar.val <- cells[[i]]
    if (bar.nreal < 1) bar.val <- 0
    bar.inf <- is.infinite(bar.val)
    bar.inf.pos <- bar.inf & bar.val > 0
    bar.inf.neg <- bar.inf & bar.val < 0
    bar.i <- sapply(bar.val, function(bi) which.min(abs(cellcolmap - bi)))
    bar.col <- names(cellcolmap)[bar.i]
    #bar.shade <- ifelse(bar.inf, 10L, NA_integer_)
    
    #bar.cell.lwd <- ifelse( bar.n < cell.limit, cell.lwd, NA)
    bar.cell.lwd = cell.lwd
    bar.shade <- ifelse(bar.inf, inf.shading, 0)
    bar.col.inf = bar.col
    if (bar.n < cell.limit) { bar.col.inf = "black" } else { bar.col.inf[bar.inf] = "black" }
    bar.col[bar.inf.neg] <- cell.col.inf[1]
    bar.col[bar.inf.pos] <- cell.col.inf[2]
    
    # omit cell borders if the bar has more than cell.limit cells
    # bar.border <- if (bar.n < cell.limit) rep("black",bar.n) else bar.col
    # make little boxes
    rect(xsteps[-(bar.n+1)], ysteps[i+1]+ygap-yspace, xsteps[-1], ysteps[i+1]+yspace,
         col = bar.col, lwd = bar.cell.lwd, border = NA) # works only whith lwd >0
    rect(xsteps[-(bar.n+1)], ysteps[i+1]+ygap-yspace, xsteps[-1], ysteps[i+1]+yspace,
         col = bar.col.inf, lwd = bar.cell.lwd, density = bar.shade) # works only whith lwd >0
    # and another box around the whole bar
    rect( xbound[1], ysteps[i+1]+ygap-yspace, xsteps[length(xsteps)], ysteps[i+1]+yspace, lwd=cell.outer )
  }
  
  
  
  title( ... )
  # XAXIS DESIGNATION
  text( (xbound[1]+xbound[2])/2, ybound[1]+0.015+strheight("0",cex = xlab.cex)*2, labels=xlab, pos=3, cex=xdes.cex )
  segments( xbound[1], ybound[1]+0.015, xbound[1], ybound[2]-yspace, lwd=cell.outer)
  
  # color legend
  if (key) {
    lc = c( 0.8,ybound[2]-ygap-yspace,0,ybound[2]-yspace*3 )
    absmax = max(abs(cellbound))
    lc.min <- min(cellbound)
    lc.max <- max(cellbound)
    cellcolmap.center <- cellcolmap[51]
    cellcolmap <- cellcolmap[lc.min <= cellcolmap & cellcolmap <= lc.max]
    lc.xsteps = seq( xbound[1], xbound[2], length.out=key.n+1 )
    lc.xgap = lc.xsteps[1] - lc.xsteps[2]
    
    lc.density <- rep(0, key.n)
    if (any(cellinf)) {
      lc.range <- c(-Inf, seq( lc.min, lc.max, length.out=key.n-2), Inf)
      lc.col <- c(cell.col.inf[1], names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n-2)], cell.col.inf[2])
      lc.density[c(1,length(lc.density))] = inf.shading
      if (key.n %% 2 > 0 && lc.min < 0 && lc.max > 0) {
        lc.range = c( -Inf, seq( lc.min, 0, length.out=(key.n-1)/2 ), seq( 0, lc.max, length.out=(key.n-1)/2 )[-1], Inf )
        print(lc.range)
      }
    } else {
      lc.range = seq( lc.min, lc.max, length.out=key.n )
      lc.col <- names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n)]
      if (key.n %% 2 > 0 && lc.min < 0 && lc.max > 0) {
        lc.range = c( seq( lc.min, 0, length.out=(key.n-1)/2+1 ), seq( 0, lc.max, length.out=(key.n-1)/2+1 )[-1] )
        lc.col <- names(cellcolmap)[seq(1,length(cellcolmap),length.out=key.n)]
      }
    }
    
    key.num = as.character(round(lc.range,1))
    if (!is.null(cell.bounds)) {
      or = range(celldata[!cellinf & !cellmis])
      if (cell.bounds[1] > min(or)) { key.num[1+(any(cellinf))] = paste0("<",key.num[1+(any(cellinf))])}
      if (cell.bounds[2] < max(or)) { key.num[length(key.num)-(any(cellinf))] = paste0(key.num[length(key.num)-(any(cellinf))],">")}
    }
    
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col=lc.col, lwd=cell.outer )
    rect( lc.xsteps[-(key.n+1)]-lc.xgap*.1, lc[2], lc.xsteps[-1]+lc.xgap*.1, lc[4], col="black", lwd=cell.lwd, density = lc.density, border=NA )
    text( (lc.xsteps[-(key.n+1)]+lc.xsteps[-1])/2, lc[2], pos=1, labels=key.num, cex=xlab.cex, font=2 )
    text( (xbound[1]+xbound[2])/2, lc[2]-strheight("0",cex=xlab.cex)*1.5 , labels=key.lab, pos=1, cex=xdes.cex )
  }
  par(xpd=T)
}
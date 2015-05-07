#' Cell Plot
#' 
#' Plots a horizontal barchart of strictly positive values in x. For each entry, a vector of additional values needs to be provided in a list. 
#' The additional values are plotted as cells subdividing the length of the corresponding bar. A typical usage scenario is plotting the 
#' enrichment of Gene Ontology terms, with individual cells reflecting the differential regulation of the constituent genes.
#' 
#' Details
#' 
#' @param x Strictly positive numeric vector.
#' 
#' @param cells List of vectors. Must be the same length as x.
#'
#' @param lab.col Vector of color names. Must be the same length as x. Defaults to black.
#' 
#' @param cell.colorFunction Function to generate the color spectrum for the cell mappings. Defaults to colorRampPalette( c("blue","white","red") ).
#'
#' @param space Scaling factor for spacing between bars.
#' 
#' @param x.mar Numeric vector defining the inside margin to the left and right of the main plot. Defaults to c(0.2,0.1).
#'
#' @param y.mar Numeric vector defining the inside margin to the top and bottom of the main plot. Defaults to c(0.08,0).
#'
#' @param lab.cex Scaling factor for label text size.
#'
#' @param xdes.cex Scaling factor for key legend text size.
#'
#' @param xlab.cex Scaling factor for x-label text size.
#'
#' @param xlab.tick Spacing between consecutive tickmarks, in the same units as used for x.
#' Defaults to max(x)/10.
#'
#' @param cell.lwd Size of the border of individual cells. Defaults to 1.
#'
#' @param cell.outer Size of outer border around bars. Defaults to 2.
#'
#' @param cell.sort Should cell values be sorted? Defaults to TRUE.
#'
#' @param cell.limit Number of cells above which separators between cells are omitted. Defaults to 50.
#'
#' @param xlab Label for x-axis.
#'
#' @param key.lab Label for color key.
#'
#' @param spacers Numeric vector. Inserts empty space after the specified positions to visually group bars.
#'
#' @param scaleTo Scale the plotting area to accomodate exactly this number of bars. Used to ensure visual consistency across different plots. 
#' Defaults to NULL, i.e. the plotting area is scaled to fit.
#'
#' @param \dots Arguments passed through to the title() function. E.g. main = "Plot title".
#' 
#' @author Robert Sehlke
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
#' 
#' ## Plot with spacers
#' cell.plot(x, cells, xcolor, spacers = c(4,8), xlab.tick = 0.5, main="Cell Plot Demo", xlab="log(enrichment)", cell.limit = 80)
#' }
#' 


#' @rdname cell.plot
#' @export
cell.plot = function( x, cells, lab.col=NULL, cell.colorFunction=colorRampPalette( c("blue","white","red") ), space=0.1, x.mar=c(0.2,0.1), y.mar = c(0.08,0), lab.cex = 1, xdes.cex=1, xlab.cex=1, xlab.tick=NULL, 
                      xlab.yoffset = 0.08, sym=T, cell.lwd=1, cell.outer=2, cell.sort=T, cell.limit=50, xlab="", key=T, key.lab="Color Key", spacers=NULL, scaleTo=NULL, ... ) {
  par(xpd=NA)
  if (is.null(xlab.tick)) { xlab.tick = round( max(x) / 10, digits = 1) }
  ticksize = xlab.tick
  ybound = c(1,0) + c(-1,1)*y.mar
  
  # scale bar area height to uniform scaleTo bars -- if it isn't provided, scale to fit
  if ( is.null(scaleTo) ) { scaleTo = length(x) + length(spacers) }
  ybound[2] = ybound[1] - ( (ybound[1] - ybound[2]) / scaleTo ) * length(x)    
  
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
  cellbound = c( min(celldata), max(celldata) )
  if (sym) { cellbound = rep( max(abs(cellbound)),2 ) * c(-1,1) }
  cellcolmap = seq( cellbound[1], cellbound[2], length.out=101 )
  names(cellcolmap) = cell.colorFunction(101)
  
  labvec = rep("",length(x))
  if (!is.null(names(x))) { labvec=names(x) }
  colvec = lab.col
  if( is.null(colvec) ) { colvec = rep("black",length(x)) }
  
  # do the actual plotting
  plot.new()
  for (i in 1:length(x) ) {
    xsteps = seq( xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1], length.out=(length(cells[[i]])+1) )
    if (is.null(cells)) { xsteps = c(xbound[1], (x[i]/max(x))*(xbound[2]-xbound[1])+xbound[1]) }
    
    text( xbound[1], ysteps[i+1]+ygap*0.5, labels=labvec[i], pos=2, cex=lab.cex, col = colvec[i] )        # row labels
    
    if(cell.sort) { cells[[i]] = sort(cells[[i]]) }
    this.n = length( cells[[i]] )
    
    text( xsteps[length(xsteps)], ysteps[i+1]+ygap*0.5, pos=4, cex=lab.cex, labels=length(cells[[i]]) )   # number of cells labels
    
    for ( j in 1:length(cells[[i]]) ) {
      this.col = names(cellcolmap)[ which(cellcolmap >= cells[[i]][j])[1] ]     # map colors to cell values
      this.cell.lwd = ifelse( length(cells[[i]]) < cell.limit, cell.lwd, 0)     # omit cell borders if the bar has more than cell.limit cells
      this.border = ifelse( length(cells[[i]]) < cell.limit, "black", this.col)
      rect( xsteps[j], ysteps[i+1]+ygap-yspace, xsteps[j+1], ysteps[i+1]+yspace, col=this.col, lwd=this.cell.lwd, border=this.border )  # make little boxes
    }
    
    rect( xbound[1], ysteps[i+1]+ygap-yspace, xsteps[length(xsteps)], ysteps[i+1]+yspace, lwd=cell.outer ) # and another box around the whole bar
  }
  
  conversion = (xbound[2]-xbound[1])/(max(x)-min(0,min(x)) )
  axis.at = seq(xbound[1], xbound[2], by=conversion*ticksize )
  axis.lab = round( seq( min(0, min(x)), (length(axis.at)-1)*ticksize, length.out=length(axis.at) ), digits=1 )
  
  axis(3, pos=ybound[1]+0.015, at = axis.at, labels = axis.lab, cex.axis=xlab.cex, padj=0.5, lwd=cell.outer )
  title( ... )
  text( (xbound[1]+xbound[2])/2, ybound[1]+0.015+strheight("0",cex = xlab.cex)*2, labels=xlab, pos=3, cex=xdes.cex )      # XAXIS DESIGNATION
  segments( xbound[1], ybound[1]+0.015, xbound[1], ybound[2], lwd=cell.outer)
  
  
  # color legend
  #lc = c( 0.8,ybound[2]-0.02-yspace*2,0,ybound[2]-yspace*2 )
  if (key) {
    lc = c( 0.8,ybound[2]-ygap-yspace*2,0,ybound[2]-yspace*2 )
    
    absmax = round( max(abs(celldata)), digits=2)
    lc.xsteps = seq( xbound[1], xbound[2], length.out=12)
    lc.xgap = lc.xsteps[1] - lc.xsteps[2]
    lc.range = round( seq( -absmax, absmax, length.out=11 ), digits=1)
    for (k in 1:length(lc.range)) {
      lc.col = names(cellcolmap)[ which(cellcolmap >= lc.range[k])[1] ]
      if (is.na(lc.col)) { lc.col = names(cellcolmap[length(cellcolmap)]) }
      rect( lc.xsteps[k]-lc.xgap*0.1, lc[2], lc.xsteps[k+1]+lc.xgap*0.1, lc[4], col=lc.col, lwd=cell.outer )  # COLOR LEGEND BOXES
      text( (lc.xsteps[k]+lc.xsteps[k+1])/2, lc[2], pos=1, labels=lc.range[k], cex=xlab.cex, font=2 )         # COLOR LEGEND TEXT
    }
    text( (xbound[1]+xbound[2])/2, lc[2]-strheight("0",cex=xlab.cex)*1.5 , labels=key.lab, pos=1, cex=xdes.cex )  # COLOR LEGEND DESIGNATION
  }
  par(xpd=T)
}


# I am a change
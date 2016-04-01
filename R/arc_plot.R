#' @title Arc Plot
#'
#' @description 
#' Plots a combined barchart for a differential expression experiment GO analysis.
#' Shows enrichment score (or similar), proportions of up- and downregulated set members, as well as
#' set overlap. Bar heights correspond to the number of genes of interest in the term. Overlap between
#' terms is visualized by an arc connecting them. The line width of the arc corresponds to the number
#' of overlapping genes, on the same scale as bar height. Arc colour scale maps to the Jaccard Index 
#' ( |Overlap| / |Union| ), bounded between 0 and 1, where black represents maximum overlap, i.e. Jaccard
#' Index of close or equal to 1. A threshold on the Jaccard index may be set to suppress arcs below the
#' specified value. Vertical ordering of terms is the result of hierarchical clustering, using hclust on
#' the pairwise Jaccard Distance matrix (1 - Jaccard Index).
#' 
#'
#' @param x Strictly positive numeric vector.
#' @param up.list List with named vectors, according to x, but only with positive 
#' (i.e. fold changes of gene expression) values.
#' @param down.list As up.list, but with negatives.
#' @param y.mar Plot area horizontal margins.
#' @param x.mar Plot area vertical margins.
#' @param x.bound tba
#' @param x.ticks tba
#' @param x.scale tba
#' @param x.gap gap between bar sections
#' @param x.arc.gap gap between rectangles and arcs
#' @param space tba
#' @param fixed.scale tba
#' @param t Threshold for Jaccard distance to display.
#' @param tc NA
#' @param main Plot title.
#' @param barcol Color of ratio bars.
#' @param cex.main Character expansion of the plot title.
#' 
#' @author 
#' Robert Sehlke [aut]\cr
#' Sven E. Templer [ctb]\cr
#' Authors of plotrix package for the plot functionality of arcs [ctb]
#'
#' @examples
#' \dontrun{
#' 
#' ### golub.deg data example
#' x = subset(golubstat, p<=.05 & significant>4 & !duplicated(genes))
#' x = head(x, 10)
#'
#' lfcvec = x$deg.log2fc
#' genevec = x$genes
#'
#' uplist = list()
#' downlist = list()
#' for (i in 1:length(lfcvec)) {
#'   names(lfcvec[[i]]) = genevec[[i]] 
#'   uplist[[i]] = lfcvec[[i]][lfcvec[[i]] > 0]
#'   downlist[[i]] = lfcvec[[i]][lfcvec[[i]] < 0]
#' }
#'
#' arc.plot( x = setNames(x$loge, x$term), up.list = uplist, down.list = downlist)
#' }
#' 

#' @export
arc.plot = function( x, up.list, down.list, y.mar=c(0,0), x.mar=c(0.5,0), x.bound=NULL, x.ticks=3, x.scale=1, x.gap = .05, x.arc.gap = .3,
                     space=0.1, fixed.scale=NULL, t=0.2, tc="black", main="GO Term Analysis",
                     barcol = c("deepskyblue2","coral"),
                     cex.main = 1.3) {
  # code for preciseArc subroutine adapted from plotrix draw.circle function (https://cran.r-project.org/web/packages/plotrix/index.html)
  preciseArc = function (x, y, radius, nv = 100, border = NULL, col = "red", lty = 1, 
                         lwd = 1e-6, startdegree=0, stopdegree=360, width = 0.1, ...) 
  {
    xylim = par("usr")
    xyasp = par("pin")
    xycr = diff(par("usr"))[c(1, 3)]
    yscale = xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
    xscale = xyasp[2]/xyasp[1] * xycr[1]/xycr[2]
    
    drange = (stopdegree - startdegree) / 360
    nv = round(nv*abs(drange))
    angle.inc = (2 * drange * pi) / nv
    angles = seq(2* pi * (startdegree/360), 2 * pi * (stopdegree/360), by = angle.inc)
    if (length(radius) < length(x)) 
      radius = rep(radius, length.out = length(x))
    if (length(col) < length(radius)) 
      col = rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
      xv1 = cos(angles) * (radius[circle] - width/2) * xscale + x[circle]
      yv1 = sin(angles) * (radius[circle] - width/2) + y[circle]
      xv2 = cos(angles) * (radius[circle] + width/2) * xscale + x[circle]
      yv2 = sin(angles) * (radius[circle] + width/2) + y[circle]
      
      polygon(c(xv1,xv2[length(xv2):1]), c(yv1,yv2[length(yv2):1]), border = NA, col = col[circle], lty = lty, 
              lwd = lwd, ...)
    }
  }
  
  # Assemble!
  up.n = sapply(up.list, length)
  down.n = sapply(down.list, length)
  all.list = list()
  props = vector(length = length(up.list))
  for (l in 1:length(up.list)) { 
    all.list[[l]] = c(up.list[[l]], down.list[[l]])
    props[l] = min( length(down.list[[l]])/length(all.list[[l]]), 1)
  }
  enrich = x
  eb = ifelse(is.null(x.bound),max(x),x.bound)
  
  # Overlap! Jaccard Indices!
  dm = matrix(0,nrow=length(all.list),ncol=length(all.list))
  jc = dm
  for (i in 1:nrow(dm)) { for (j in 1:ncol(dm)) {
    t1 = all.list[[i]]
    t2 = all.list[[j]]
    dm[i,j] = length(intersect(t1, t2))
    jc[i,j] = dm[i,j] / length(union(t1,t2))
  }}
  dimnames(dm) = list(names(enrich),names(enrich))
  
  # Cluster! (using the inverted Jaccard similarity)
  cl = hclust(as.dist(1-jc))$order
  overlapmatrix = dm[cl,cl]
  jcmatrix = jc[cl,cl]
  props=props[cl]
  enrich=enrich[cl]
  
  # Aspect ratios! (give partial credit to Plotrix developers)
  xylim = par("usr")
  xyasp = par("pin")
  xycr = diff(par("usr"))[c(1, 3)]
  yscale = xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
  xscale = xyasp[2]/xyasp[1] * xycr[1]/xycr[2]
  yscale = (diff(par("usr")[3:4])/par("pin")[2])
  
  # Define plotting area and scaling and stuff!
  par(xpd=NA)
  ybound = c(1,0) + c(-1,1)*y.mar
  xbound = c(0,1) + c(1,-1)*x.mar
  if ( !is.null(fixed.scale) ) {
    ybound[2] = ybound[1]-( 0.3 * yscale * fixed.scale * nrow(overlapmatrix) )
    ysteps = ybound[2] + c(0, cumsum( rep(0.3 * yscale * fixed.scale, nrow(overlapmatrix) ) ) )
    
    if ( ybound[2] < par("usr")[3] ) {
      warning("Plotting area too small! Decrease fixed.scale or increase space.")
    }
  }
  ysteps = seq( ybound[2], ybound[1], length.out=( nrow(overlapmatrix)+1 ) )
  ygap = abs(ysteps[1]-ysteps[2])
  yspace = space * ygap * 0.5
  
  # Scale and threshold overlap matrix!
  overlapmatrix = overlapmatrix/max(overlapmatrix)
  overlapmatrix = overlapmatrix * (jcmatrix > t) * ygap * (1-space)
  rldiag = diag(overlapmatrix)
  
  # arc color matrix
  acn = 100
  arcol = seq(0,1,length.out = 10)
  names(arcol) = colorRampPalette(c("white","black"))(10)
  arcolmat = apply( jcmatrix, c(1,2), function(x) x = names(arcol)[which.min(x >= arcol)] )
  
  
  # PLOTTING OVER HERE!
  
  gapex = 4 * x.scale
  ticks = x.ticks
  plot.new()
  
  # arcs
  ao = order(jcmatrix)
  aov = jcmatrix[ao]
  for (o in 1:length(aov) ) {
    ti = which(jcmatrix == aov[o], arr.ind = T)
    j = ti[1,1]
    k = ti[1,2]
    if ( (overlapmatrix[j,k] > 0) && (j != k)) {
      p1 = mean(ysteps[j:(j+1)])
      p2 = mean(ysteps[k:(k+1)])
      hinge = mean( c(p1,p2) )
      rad = abs(p1-p2)/2
      preciseArc(x.mar[1]+x.arc.gap*ygap*gapex*xscale*0.1,hinge,rad,width = overlapmatrix[j,k], col = arcolmat[j,k], startdegree = -90, stopdegree = 90, nv = 100)
    }
  }
  
  
  # enrichment axis ticks
  at.enrich = seq( x.mar[1]-ygap*gapex*xscale*2,
                   x.mar[1]-ygap*gapex*xscale*(1+x.gap),
                   length.out = ticks)
  lab.enrich = round( seq( eb, 0, length.out = ticks), digits = 1 )
  egap = diff(at.enrich[c(1,length(at.enrich))])
  
  # up / down bars
  rect( x.mar[1]-ygap*gapex*xscale , 
        ysteps[1:(length(ysteps)-1)]+(ygap-rldiag)/2, 
        x.mar[1]-ygap*gapex*xscale*(1-props), 
        ysteps[2:(length(ysteps))]-(ygap-rldiag)/2, 
        col=barcol[1], border=NA)
  rect( x.mar[1]-ygap*gapex*xscale*(1-props) ,
        ysteps[1:(length(ysteps)-1)]+(ygap-rldiag)/2, 
        x.mar[1], 
        ysteps[2:(length(ysteps))]-(ygap-rldiag)/2, 
        col=barcol[2], border=NA)
  
  # enrichment grid and bars
  segments(at.enrich, ybound[2]-ygap*space, at.enrich, ybound[1], col = "darkgrey")
  segments(at.enrich[length(at.enrich)], ybound[1], at.enrich[1], ybound[1], col = "darkgrey")
  rect( at.enrich[1] + egap*(1-enrich/eb),
        ysteps[1:(length(ysteps)-1)]+(ygap-rldiag)/2,
        at.enrich[length(at.enrich)], 
        ysteps[2:(length(ysteps))]-(ygap-rldiag)/2, 
        col="black", border = NA)
  
  
  # axis with custom axis labels
  axis(1, pos=ybound[2]-ygap*space, 
       at = at.enrich, 
       labels=NA, 
       cex.axis=ifelse(is.null(fixed.scale),1,fixed.scale), lwd = 1, padj = ifelse(is.null(fixed.scale), -1, -1/fixed.scale*yscale ))
  text( at.enrich, 
        ybound[2]-ygap*1.5-ygap*space, 
        lab.enrich, 
        cex=ifelse(is.null(fixed.scale),1,fixed.scale))
  
  
  
  # color legend
  rect( x.mar[1]-ygap*gapex*xscale , ybound[2]-ygap, 
        x.mar[1]-ygap*gapex*xscale*0.75, ybound[2]-ygap*space, 
        col=barcol[1], border = NA)
  rect( x.mar[1]-ygap*gapex*xscale , ybound[2]-ygap*2, 
        x.mar[1]-ygap*gapex*xscale*0.75, ybound[2]-ygap*space-ygap, 
        col=barcol[2], border = NA)
  
  # color legend labels
  text( x.mar[1]-ygap*gapex*xscale*0.6+strwidth(c("% Downregulated","% Upregulated"),cex=ifelse(is.null(fixed.scale),1,fixed.scale)/2), 
        c( ybound[2]-ygap*0.5-ygap*space, ybound[2]-ygap*1.5-ygap*space), 
        c("% Downregulated","% Upregulated"), 
        cex=ifelse(is.null(fixed.scale),1,fixed.scale))
  
  # axis label
  text( mean(at.enrich),
        ybound[2]-ygap*2.5-ygap*space,
        "Enrichment",
        cex=ifelse(is.null(fixed.scale),1,fixed.scale) )
  
  # row labels
  text( x.mar[1]-ygap*gapex*xscale*2, 
        (ysteps[1:(length(ysteps)-1)] + ysteps[2:(length(ysteps))])/2, 
        rownames(overlapmatrix), pos = 2, 
        cex=ifelse(is.null(fixed.scale),1,fixed.scale) )
  
  # main title
  text( x.mar[1]-ygap*gapex*xscale*2, 
        ybound[1]+ygap*xscale*2, 
        main, 
        cex=ifelse(is.null(fixed.scale),cex.main,fixed.scale*cex.main) )
  # return(overlapmatrix)
}
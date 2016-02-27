#' @title GO Histogram
#'
#' @description 
#' Plots a histogram, combining multiple topGO analyses into one figure. Bars
#' show the absolute numbers of up- and downregulated genes. Terms that were
#' detected as significantly enriched in the respective analysis are coloured,
#' grey bars indicate no enrichment. In addition, p-values from a one-sided
#' t-test on the significant log fold changes may be included, and are visualized
#' by asterisks at the end of the bars.
#'
#' @param framelist List of TopGO data frames. The frames have to include two
#' additional columns: 'Upregulated' and 'Downregulated', with the numbers of
#' up- and down-regulated genes, respectively.
#' 
#' @param alpha Significance threshold for enrichment. Terms that are not
#' significantly enriched are still displayed using grey bars. Defaults to 0.05.
#' 
#' @param alpha.term Identifies which column contains the p-values of the enrichment
#' analysis. Defaults to "Elim".
#' 
#' @param min.sig Minimum times a term has to be detected as significantly
#' enriched across all groups to be displayed in the histogram. Defaults to 1.
#'
#' @param min.genes Minimum number of genes for a significantly enriched term
#' to be included. Additional to the min.sig filter. Defaults to 10.
#' 
#' @param max.genes Maximum number of genes for a significantly enriched term
#' to be included. Defaults to 100.
#' 
#' @param bar.scale If not NULL (the default), set bar height to a fixed value relative to the
#' numeric factor provided. Use to ensure consistency across plots with varying numbers of elements.
#' Returns a warning if the plotting area is too small and uses the relative scale instead.
#' 
#' @param reorder If TRUE, terms significantly enriched in the same groups appear together.
#' Defaults to TRUE.
#' 
#' @param main Plot title.
#' 
#' @param show.go.id If TRUE, shows the GO ID number along with the term.
#' 
#' @param show.ttest If TRUE, displays p-values from one-sided t-tests beside the bars, coded as
#' asterisks. The p-values need to be included in the input data frames as columns "p.up.adj" 
#' and "p.down.adj" (NOTE: make more generic!). Defaults to FALSE.
#' 
#' @param lab.cex Scalar, Character expansion factor for labels.
#' 
#' @param axis.cex Scalar, Character expansion factor for axis tick labels.
#' 
#' @param group.cex Scalar, Character expansion factor for group labels.
#' 
#' @param go.selection Character vector of GO IDs that should be displayed.
#' 
#' @param term.selection Character vector of GO term identifiers that should be displayed.
#'
#' @author 
#' Robert Sehlke [aut]\cr
#' Sven E. Templer [ctb]
#'
#' @examples
#' \dontrun{
#' 
#' }
#'

#' @export
go.histogram = function( framelist, alpha=0.05, alpha.term="Elim", min.sig=1, min.genes=10, max.genes=100, bar.scale=NULL,
                         reorder=T, main="GO enrichment", show.go.id=FALSE, prefix="",show.ttest=F, lab.cex=1, 
                         axis.cex=1, group.cex=NULL, go.selection=NULL, term.selection=NULL) {
  
  # Auxiliary functions -- put elsewhere?
    addsig = function( f, sigcolumn, at ) {
      tt = rep("",nrow(f))
      tt[which(f[,sigcolumn] < 0.05)] = "*"
      tt[which(f[,sigcolumn] < 0.01)] = "**"
      tt[which(f[,sigcolumn] < 0.001)] = "***"
      ttsw = strwidth(tt)/1.5
      text(f[,grep(paste0(prefix,"Upregulated"),colnames(f))] + ttsw, at, tt)
    }
    
    grep.multi = function(xvec, target) {
      acc = vector(mode="numeric")
      for (i in 1:length(xvec)) {
        acc = c(acc, grep(xvec[i], target))
      }
      return(unique(acc))
    }
  
  # Axis limit (for now, symmetric)
    lim = max( sapply(framelist, function(x) max(x[,paste0(prefix, c("Upregulated","Downregulated"))])  ) )
  
  # preprocess into one combined data frame
    if (!is.null(go.selection)) {
      for (f in 1:length(framelist)) {
        framelist[[f]] = framelist[[f]][which(framelist[[f]]$GO.ID %in% go.selection),]
      }
    }
    if (!is.null(term.selection)) {
      for (f in 1:length(framelist)) {
        framelist[[f]] = framelist[[f]][which(framelist[[f]]$Term %in% term.selection),]
      }
    }
    if(is.null(names(framelist))) { names(framelist) = paste0("group_",1:length(framelist) ) }
    framelist = lapply(framelist, function(x) { x$idterm = paste(x$GO.ID, x$Term); return(x) })
  
  
  # Integrate the topGO frames into one combined data.frame
    allrows = vector()
    for (i in 1:length(framelist)) { 
      sel = which(!is.na(framelist[[i]][,"idterm"]))
      framelist[[i]] = framelist[[i]][sel,]
      rownames(framelist[[i]]) = framelist[[i]][,"idterm"]
      allrows = c(allrows, rownames(framelist[[i]]))
    }
    allrows = unique(allrows)
    
    intframe = data.frame(row.names=allrows)
    outcolnames = vector()
    for (j in 1:length(framelist)) {
      tmp = data.frame( framelist[[j]] )
      colnames(tmp) = colnames( framelist[[j]] )
      rownames(tmp) = rownames(framelist[[j]])
      outcolnames = c( outcolnames, paste( names(framelist[j]), colnames(tmp), sep="." ) )
      intframe = data.frame( cbind(intframe, tmp[rownames(intframe),]) )
      rownames(intframe) = allrows
    }
    colnames(intframe) = outcolnames
  
  
  # select terms if specified
    if (!is.null(term.selection)) {
      termorder = gsub( "GO:[0-9]+ ","", rownames(intframe) )
      termorder = match(term.selection, termorder)
      intframe = intframe[rev(termorder),]
    }
  
  # apply filters
    filterframe = (intframe[,grep(paste0("\\.", alpha.term),colnames(intframe))] <= alpha) *  ( intframe[,grep("\\.Significant",colnames(intframe))] >= min.genes )
    intframe = intframe[ apply(filterframe,1,function(x) sum(x, na.rm = T)) >= min.sig, ]
    intframe = intframe[ which(apply( intframe[,grep("\\.Significant",colnames(intframe))], 1, function(x) all( x <= max.genes, na.rm = T) ) ), ]
    # select terms with at least min.sig samples
    
  # cluster by significance status (p[enrich]<=alpha) across groups (other options?)
    if (reorder) { 
      sf = intframe[,grep(paste0("\\.", alpha.term),colnames(intframe))]
      sf = sf <= alpha
      sf[is.na(sf)] = FALSE
      sfo = hclust(dist(sf))$order
      intframe = intframe[sfo,] 
    }
  
  # limit for axis
    lim = max(intframe[,grep.multi(paste0(prefix, c("Upregulated","Downregulated")),colnames(intframe))], na.rm=T)
  
  # Layout
    dv = dev.size(units = "in")
    m = m.old = par("mai")
    m[c(1,3)] = sapply( m[c(1,3)]/dv[2], function(x) max(x, 0.01) )
    m[c(2,4)] = sapply( m[c(2,4)]/dv[1], function(x) max(x, 0.01) )
    nc = length(framelist) * 2      # n.o. columns, just a shortcut variable
    thi = 0.6; th = thi / dv[2]                # desired label area height in inches (then scaled)
    bhi = 0.18; bh = (bhi * ifelse(is.null(bar.scale),1,bar.scale))  / dv[2]       # height of one bar if fixed in inches
    ta = 2 / dv[1]    # label area in inches (then scaled)
    bot = 0           # bottom margin of bar plots
    top = 0.4         # top margin of bar plots (need to accomodate the axes)
    gap = 0.15         # space to the left/right of bars
  
    if(!is.null(bar.scale)) { 
      neededspace = bhi * nrow(intframe) + sum( m.old[c(1,3)] ) + thi
      if ( neededspace > dv[2] ) { stop(paste0("figure region too small (set device height to above ",round(neededspace,digits = 2)," in)")) }
    }
  
    layoutmatrix = rbind( c( nc/2+nc+5, rep( nc/2+nc+5, nc+1) ),                # m[3] # top margin  
                          c( nc/2+nc+4, sort( rep( 1:(nc/2), 2) ), nc/2+1 ),    # th
                          c( nc/2+nc+4, (nc/2+2):(nc/2+nc+1), nc/2+nc+2 ),      # bararea
                          c( nc/2+nc+3, rep( nc/2+nc+3, nc+1) ) )               # m[1] # bottom margin
  
    barareaheight = ifelse( is.null(bar.scale), 
                            1-m[3]-m[1]-th, 
                            min( (top+bot)/dv[2]+bh*nrow(intframe), 1-m[3]-m[1]-th ) )
    barplotwidth = (1-m[2]-m[4]-ta)/nc
  
    if ( any( c(barareaheight, barplotwidth-gap/dv[1])  <= 0 ) ) { par(mai=m.old, xpd = F); layout(1); stop("figure region too small") }
  
  layout(layoutmatrix, 
           heights = c( m[3],th,barareaheight,1-m[3]-th-barareaheight),
           widths = c( m[2], rep( barplotwidth, nc ), m[4]+ta ) )
  
  print(m)
  print(c( m[3],th,barareaheight,1-m[3]-th-barareaheight))
  print(c( m[2], rep( barplotwidth, nc ), m[4]+ta ))
  
  # Group labels and title
    par(mai=c(0,0,0,0), xpd = T)
    for (i in 1:length(framelist) ) {
      plot.new()
      text(0.5,0.5,names(framelist)[i], cex = ifelse(is.null(group.cex), 1.5, group.cex) )
    }
    plot.new()
    text(0.5,0.5,main, cex=1)
  
  # Cycling through the groups and plotting each in turn
  for (i in 1:length(framelist) ) {
    tmp = intframe[,grep(names(framelist[i]), colnames(intframe))]
    ts = tmp[,grep(paste0("\\.", alpha.term),colnames(tmp))]
    ucol = rep("darkgrey",length(ts))
    dcol = rep("grey",length(ts))
    ucol[which(ts <= alpha)] = "coral"
    dcol[which(ts <= alpha)] = "deepskyblue2"
    
    
    # left-pointed bars / downregulated
      par(mai=c(bot,gap,top,0))
      
      vline = barplot(-tmp[,grep(paste0(prefix,"Downregulated"),colnames(tmp))], horiz = T, 
                      xlim = c(-lim,0), border=NA, col = dcol, axes=F)       # getting the bar positions
      axis(3, at=-round( seq(0,lim,length.out = 3)), labels = round( seq(0,lim,length.out = 3)), cex.axis=axis.cex, las=3 )
      abline(h=vline, col="grey", lty=2)
      # plotting the barplot a second time to overlay the marker lines -- any better idea how to do this?
      barplot(-tmp[,grep(paste0(prefix,"Downregulated"),colnames(tmp))], horiz = T, 
              xlim = c(-lim,0), border=NA, col = dcol, axes=F,add=T)
      
      # adding significance asterisks
      if(show.ttest) { addsig( tmp, paste0(names(framelist)[i],".p.down.adj"), vline) }
    
    
    # right-pointed bars / upregulated, and middle lines
      par(mai=c(bot,0,top,gap))
      barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
              xlim = c(0,lim), border=NA, col=ucol, axes=F)
      abline(h=vline, col="grey", lty=2)
      barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
              xlim = c(0,lim), border=NA, col=ucol, axes=F, add=T)
      axis(3, at=round( seq(0,lim,length.out = 3)), cex.axis=axis.cex, las=3 )
      lines(c(0,0),c( par('usr')[4], -0.5*vline[1] ) )
      
      if(show.ttest) { addsig( tmp, paste0(names(framelist)[i],".p.up.adj"), vline) } 
    
  }
  
  # Row labels
    par(mai=c(bot,0,top,0))
    barplot( rep(NA,length(vline)), axes = F, horiz=T, xlim = c(0.1,1) ) # dummy plot to get the scale -- any better way?
    text(0.1, vline, pos = 4, cex=lab.cex,
         labels = if (show.go.id) { rownames(intframe) } else { gsub("GO:[0-9]+ ","",rownames(intframe)) })

  # bottom margin and left margin placeholder plots
    par(mar=c(0,0,0,0) )
    plot.new()
    plot.new()
    plot.new()
  
  # Reset graphical params
    par(mai=m.old, xpd = F)
    layout(1)
  
}
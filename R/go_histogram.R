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
go_histogram = function( framelist, alpha=0.05, min.sig=1, min.genes=10, max.genes=100, scaleTo=NULL, reorder=T,
                                   golabels="idterm", main="GO enrichment", show.go.id=FALSE, prefix="",show.ttest=F, lab.cex=1, axis.cex=1,
                                   go.selection=NULL, term.selection=NULL, scale.sep=1) {
  
  addsig = function( f, sigcolumn, at ) {
    tt = rep("",nrow(f))
    tt[which(f[,sigcolumn] < 0.05)] = "*"
    tt[which(f[,sigcolumn] < 0.01)] = "**"
    tt[which(f[,sigcolumn] < 0.001)] = "***"
    ttsw = strwidth(tt)/1.5
    text(f[,grep(paste0(prefix,"Upregulated"),colnames(f))] + ttsw, at, tt)
  }
  
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
  
  
  
  # Integrate the topGO frames
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
  filterframe = (intframe[,grep("\\.Elim",colnames(intframe))] <= alpha) *  ( intframe[,grep("\\.Significant",colnames(intframe))] >= min.genes )
  intframe = intframe[ apply(filterframe,1,function(x) sum(x, na.rm = T)) >= min.sig, ]
  intframe = intframe[ which(apply( intframe[,grep("\\.Significant",colnames(intframe))], 1, function(x) all( x <= max.genes, na.rm = T) ) ), ]
  
  # cluster by significance across groups
  sf = intframe[,grep("\\.Elim",colnames(intframe))]
  sf = sf <= alpha
  sf[is.na(sf)] = FALSE
  sfo = hclust(dist(sf))$order
  
  if (reorder) { intframe = intframe[sfo,] }
  if( is.null(scaleTo) ) { scaleTo = nrow(intframe) }
  if( nrow(intframe) > scaleTo) { scaleTo = nrow(intframe) }
  
  print(nrow(intframe))
  #return(intframe)
  # hack to uniformly scale to a specific number of max rows
  trd = nrow(intframe)
  filler = matrix(NA, nrow= scaleTo - trd, ncol=ncol(intframe))
  used = (scaleTo - trd + 1):scaleTo
  colnames(filler) = colnames(intframe)
  intframe = rbind(filler, intframe)
  
  # limit for axis
  lim = max(intframe[,grep.multi(paste0(prefix, c("Upregulated","Downregulated")),colnames(intframe))], na.rm=T)
  
  
  # The Actual Plotting
  #plot.new()
  nc = length(framelist) * 2
  th = 1    # text height
  ta = 0.3  # text area as fraction of total width
  layoutmatrix = rbind( c( sort( rep( 1:(nc/2), 2) ), rep( nc/2+1, 1 ) ),
                        c( (nc/2+2):(nc/2+nc+1), rep(nc/2+nc+2, 1 ) ) )
  devheight = ifelse( th/dev.size(units = "cm")[2] < 1, th/dev.size(units = "cm")[2], 0.2 )
  layout(layoutmatrix, 
         heights = c(devheight,1-devheight),
         widths = c( rep( (1-ta)/nc, nc ), ta ) )
  
  
  par(mar=c(0,0,0,0))
  for (i in 1:length(framelist) ) {
    plot.new()
    text(0.5,0.5,names(framelist)[i], cex = 1.5)
  }
  plot.new()
  text(0.5,0.5,main, cex=1)
  
  par( xpd=T )
  
  for (i in 1:length(framelist) ) {
    tmp = intframe[,grep(names(framelist[i]), colnames(intframe))]
    ts = tmp[,grep("\\.Elim",colnames(tmp))]
    ucol = rep("darkgrey",length(ts))
    dcol = rep("grey",length(ts))
    ucol[which(ts <= alpha)] = "coral"
    dcol[which(ts <= alpha)] = "deepskyblue2"
    
    
    # left-pointed bars
    par(mar=c(3,1,2,0))
    
    vline = barplot(-tmp[,grep(paste0(prefix,"Downregulated"),colnames(tmp))], horiz = T, 
                    xlim = c(-lim,0), border=NA, col = dcol, axes=F)
    axis(3, at=-round( seq(0,lim,length.out = 3)), labels = round( seq(0,lim,length.out = 3)), cex.axis=axis.cex, las=3 )
    abline(h=vline[used], col="grey", lty=2)
    barplot(-tmp[,grep(paste0(prefix,"Downregulated"),colnames(tmp))], horiz = T, 
            xlim = c(-lim,0), border=NA, col = dcol, axes=F,add=T)
    
    if(show.ttest) { addsig( tmp, paste0(names(framelist)[i],".p.down.adj"), vline) }
    
    # right-pointed bars and middle lines
    par(mar=c(3,0,2,1))
    barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
            xlim = c(0,lim), border=NA, col=ucol, axes=F)
    abline(h=vline[used], col="grey", lty=2)
    barplot(tmp[,grep(paste0(prefix,"Upregulated"),colnames(tmp))], horiz = T, 
            xlim = c(0,lim), border=NA, col=ucol, axes=F, add=T)
    axis(3, at=round( seq(0,lim,length.out = 3)), cex.axis=axis.cex, las=3 )
    lines(c(0,0),c( par('usr')[4], vline[min(used)]-vline[1] ) )
    
    if(show.ttest) { addsig( tmp, paste0(names(framelist)[i],".p.up.adj"), vline) } 
    
  }
  
  # labels
  par(mar=c(3,0,2,0))
  barplot( rep(0,length(vline[used])), axes = F, horiz=T, xlim = c(0.1,1) )
  abline(h = 0.1, col="white")
  text(0.1,vline[used], , pos = 4, cex=lab.cex,
       labels = if (show.go.id) { rownames(intframe)[used] } else { gsub("GO:[0-9]+ ","",rownames(intframe)[used]) })
  
}
# go_hist2.r

gohist <- function (
  x, 
  name.label = c("GO.ID","Term"),
  name.loge  = "LogEnrich",
  name.p     = "pvalCutOff",
  name.logfc = "log2FoldChange",
  name.px    = "pvalue",
  alpha      = 0.05,
  alpha.deg  = 0.05,
  sub.alpha  = 1,
  sub.ids    = NULL,
  col.sign   = "seagreen",
  col.norm   = "grey70",
  col.updown = c("coral","deepskyblue2"),
  par.mar    = c(3,1,3,1),
  size.label = 1,
  ...
  )
{
  
  ## check
  if (!all(sapply(x, function (i) c(name.label, name.loge, name.p, name.logfc, name.px) %in% colnames(i))))
    stop("missing columns")
  par0 <- par()
  name.label <- name.label[1] # match.arg?
  if (is.null(names(x)))
    names(x) <- paste("group", seq_along(x))
  
  ## preformat
  
  x <- lapply(x, function (i) {
    i <- i[,c(name.label, name.loge, name.p, name.logfc, name.px)]
    names(i) <- c("id","e","ep","f","fp")
    return(i)
  } )
  
  ## select
  
  if (!is.null(sub.ids)) x <- lapply(x, subset, id %in% sub.ids)
  x <- lapply(x, function (i) i[i$p <= sub.alpha,])
  ids <- sort(unique(unlist(lapply(x, "[[", "id"))))
  ids <- data.frame(id = ids, stringsAsFactors = F)
  emax <- max(unlist(sapply(x, "[[", "e")), na.rm = T)
  col.updown <- colorRampPalette(col.updown)(10)
  
  ## expand
  
  x <- lapply(x, merge, x = ids, all.x = T)
  x <- lapply(x, function (i) {
    mis <- is.na(i$e)
    i$e[mis] <- 0
    i$ep[mis] <- 1
    #i$fp[mis] <- numeric()
    #i$f[mis] <- numeric()
    i$col <- ifelse(i$ep <= alpha, col.sign, col.norm)
    #i$colx <- sapply(i$f, function (j) sum(j>0)/length(j))
    #i$colx <- col.updown[as.integer(i$colx * 10)]
    i$n <- Map(function(f,p){
      u <- f>0
      d <- f<0
      su <- sum(p[u]<=alpha.deg)
      sd <- sum(p[d]<=alpha.deg)
      paste0(su, "\n", sn)
    },i$f,i$fp)
    i$n[mis] <- ""
    return(i)
  })
  ## dimensions
  
  nr <- nrow(ids)
  nc <- length(x) * 2 + 1
  m <- matrix(seq(nc), 1, nc)
  #wc <- rep(1, nc)
  wc <- c(10, rep(c(3,10), (nc-1)/2))
  #return(wc)
  wc[1] <- wc[1] * size.label
  layout(m, widths = wc)
  par(mar = par.mar)
  
  ## plot
  
  # axis
  y <- barplot(rep(1,nr), horiz = T, border = NA, axes = F, col = "grey90", main = "\nTerm")
  text(1, y, ids$id, adj = 1)
  
  # bars
  null <- Map(function (i,n) {
    y <- barplot(rep(1,nr), horiz = T, xlim = c(-.1,1.1), col = "grey80", axes = F, main = "\nn(DEG)") #i$colx
    text(.5, y, i$n)
    barplot(i$e, horiz = T, xlim = c(0, emax), col = i$col, main = n) # , names.arg = i$n, las = 2
  }, x, names(x))
  
  #par(par0)
  return(x) # invisible(NULL)
  
}
### build script for the CellPlot package

# set working directory to package source folder root

### package

library(devtools)
library(rmarkdown)
render('vignettes/CellPlotManual.Rmd')
document()
install(build_vignettes = TRUE)
library(CellPlot)
?CellPlot
vignette("CellPlotManual")

### paper image

library(CellPlot)
data(golubstat)
x <- subset(golubstat, p<=0.05 & significant>4 & !duplicated(genes))
x <- head(x, 10)
x$deg.up <- lapply(Map(setNames, x$deg.log2fc, x$genes), function (i) { i[i>0] })
x$deg.down <- lapply(Map(setNames, x$deg.log2fc, x$genes), function (i) { i[i<0] })

pdf("/beegfs/group_cd/data/projects/departments/Christoph_Dieterich/Paper_CellPlot/cpex_4_with_arc.pdf", 7, 10)

layout(matrix(1:3,nrow=3))
par(mar=c(3,0,4,0))

par(usr=c(0,1,0,1))
cell.plot( x = setNames(x$loge, x$term), cells = x$deg.log2fc, 
           y.mar = c(0, 0), x.mar = c(0.3, 0), bar.scale = .1, 
           main ="",
           x.bound = 2, cell.bounds = c(-6, 6), space = .15, key.n = 7, cex = 1.6, cell.outer = 3 )
text(0, 1.1, "A", cex=2)

par(usr=c(0,1,0,1))
sym.plot( x = setNames(x$loge, x$term), cells = x$deg.log2fc, x.annotated = x$annotated,
          y.mar = c(0.3, 0), x.mar = c(0.284, 0), 
          key.lab = "Enrichment", ticksize = 5, 
          key.n = 7, cex = 1.6, axis.cex = 1, group.cex = 1, mid.cex = 1)
text(0, 1.1, "B", cex=2)

par(mar=c(3,0,4,0))
par(usr=c(0,1,0,1))
arc.plot( x = setNames(x$loge, x$term), up.list = x$deg.up, down.list = x$deg.down, 
          y.mar = c(0, 0.1), x.mar = c(.95, 0.2), x.scale = 2.6,
          main = "", x.ticks = 5, x.bound = 2 )
text(0, 1.1, "C", cex=2)

dev.off()


### test .Rbuildignore

cat("foo\n")

#' @export
foo <- function () { NULL }


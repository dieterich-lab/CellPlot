### build script for the CellPlot package

stop("dontrun")

########################
### package
########################

rmarkdown::render('vignettes/CellPlotManual.Rmd', rmarkdown::html_document(toc = TRUE, highlight = "tango"))
#rmarkdown::render('vignettes/CellPlotManual.Rmd', rmarkdown::pdf_document(toc = TRUE, highlight = "tango"))
#rmarkdown::render('vignettes/CellPlotManual.Rmd', rmarkdown::pdf_document(toc = TRUE))
devtools::document()
devtools::install(build_vignettes = T)
library(CellPlot)
?CellPlot
vignette("CellPlotManual")

########################
### test .Rbuildignore
########################

cat("foo\n")

#' @export
foo <- function () { NULL }

########################
### dataset golubGO
########################

# CRAN
library(parallel)
library(plyr)
library(miscset)
# Bioconductor
library(topGO)
library(annotate)
library(hu6800.db)
library(multtest)

data(golub)

A <- split(golub, seq(nrow(golub)))
DEG <- do.call(rbind, mclapply(A, function(x){
  t <- t.test(x[as.logical(golub.cl)], x[!as.logical(golub.cl)])
  t$fc <- t$estimate[1]/t$estimate[2]
  suppressWarnings(t$fc <- log2(t$fc))
  data.frame(
    mean.aml = t$estimate[1],
    mean.all = t$estimate[2],
    log2fc  = t$fc,
    p = t$p.value,
    row.names = NULL)
}, mc.cores = 3))
DEG$padj <- p.adjust(DEG$p)
DEG <- data.frame(gene = golub.gnames[,3], DEG, stringsAsFactors = F)
DEG <- subset(DEG, !is.na(log2fc) & !is.na(padj))

M <- select(hu6800.db, DEG$gene, c("PROBEID","GO"))
genes <- dlply(subset(M, ONTOLOGY == "BP"), "PROBEID", function (x) unique(x$GO))
DEG <- subset(DEG, gene %in% names(genes))

## topGO object
GO <- new("topGOdata", ontology = "BP", description = 'golub',
          allGenes = setNames(DEG$padj, DEG$gene),
          geneSelectionFun = function (allScore) { allScore <= 0.05 },
          annotationFun = annFUN.gene2GO, gene2GO = genes)

GOsig <- lapply(list(golub=GO), function (x) {
  t <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  s <- getSigGroups(x, t)
  GenTable(x, pvalCutOff = s, topNodes = length(x@graph@nodes))
})

golubGO <- Map(CellPlot::mergeGOdeg, GOsig, list(DEG), list(M), map.gene = "PROBEID", deg.stats = c("padj","log2fc"))
golubGO <- lapply(golubGO, subset, !is.na(PROBEID))
save(golubGO, file = "data/golubGO.rdata")

########################
### dataset leukemiasGO
########################

# CRAN
library(parallel)
library(plyr)
library(miscset)
# Bioconductor
library(BiocParallel)
library(DESeq2)
library(topGO)
library(annotate)
library(hu6800.db)
library(leukemiasEset)

data(leukemiasEset)

M <- select(hu6800.db, featureNames(leukemiasEset), c("ENSEMBL","GO"), keytype = "ENSEMBL")
M <- subset(M, ONTOLOGY == "BP", c("ENSEMBL","GO"))
M <- M[!duplicated(M),]

genes <- dlply(M, "ENSEMBL", function (x) unique(x$GO))

A <- as(leukemiasEset,"data.frame")
A <- subset(A, LeukemiaType %in% c("NoL","ALL", "AML", "CLL"))
s <- subset(A, select = "LeukemiaType")
s$LeukemiaType <- relevel(factor(as.character(s$LeukemiaType)), "NoL")
A <- subset(A, select = colnames(A) %in% names(genes))
A <- sapply(A, as.integer)
A <- t(A)

DEG <- DESeqDataSetFromMatrix(countData = A, colData = s, design = ~ LeukemiaType)
DEG <- DESeq(DEG, fitType = "mean") #, parallel = T, BPPARAM = MulticoreParam(3))
DEG <- lapply(levels(s$LeukemiaType)[-1], function (x) {
  x <- results(DEG, c("LeukemiaType", x, "NoL"), "LeukemiaType", alpha = .05)
  x <- as.data.frame(x)
  x$gene <- rownames(x)
  return(x)
})
names(DEG) <- levels(s$LeukemiaType)[-1]

GO <- lapply(DEG, function (x) new(
  "topGOdata", ontology = "BP", description = 'Leukemia', allGenes = setNames(x$pvalue, x$gene), 
  geneSelectionFun = function (allScore) { allScore <= 0.05 }, annotationFun = annFUN.gene2GO, gene2GO = genes))

GOsig <- mclapply(GO, function (x) {
  t <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  s <- getSigGroups(x, t)
  GenTable(x, pvalCutOff = s, topNodes = length(x@graph@nodes))
}, mc.cores = 3)

leukemiasGO <- Map(CellPlot::mergeGOdeg, GOsig, DEG, list(M), map.gene = "ENSEMBL")
leukemiasGO <- lapply(leukemiasGO, subset, !is.na(ENSEMBL))
save(leukemiasGO, file = "data/leukemiasGO.rdata")

########################
### paper figure 1 (leukemiasGO version)
########################

library(CellPlot)
library(miscset)
data(leukemiasGO)

### ERRORS IN ARC-PLOT
pdf("~/tmp/fig1b.pdf", height = 10) #svg("~/tmp/fig1.svg", height = 10)
x <- leukemiasGO$CLL
x <- subset(x, Annotated > 5 & sapply(log2FoldChange, function(i)any(i>0)) & sapply(log2FoldChange, function(i)any(i<0)))
#& pvalCutOff <= 0.05 
x$Enrichment <- x$Significant / x$Expected
x$log2FoldChange <- Map(setNames, x$log2FoldChange, x$ENSEMBL)
x$Down <- lapply(x$log2FoldChange, function (y) y[y<0])
x$Up <- lapply(x$log2FoldChange, function (y) y[y>0])
x <- sort(x, TRUE, "Enrichment")
x <- head(x, 8)
layout(matrix(1:3,nrow=3))
#par(mar=c(3,0,4,0))
#par(usr=c(0,1,0,1))
cell.plot(x = setNames(x$Enrichment, x$Term), cells = x$log2FoldChange, x.mar = c(0.3,0), y.mar = c(0.2,0.1),
          space = .2, bar.scale = .7,
          main = "GO enrichment in CLL/NoL differential gene expression")
text(0, 1.1, "A", cex=2)
#par(usr=c(0,1,0,1))
sym.plot(x = setNames(x$Enrichment, x$Term), cells = x$log2FoldChange, x.annotated = x$Annotated,
         x.mar = c(0.3,0), y.mar = c(0.3, 0.0), cex = 1.6, ticksize = 5, key.lab = "Enrichment")
text(0, 1.1, "B", cex=2)
#par(mar=c(3,0,4,0))
#par(usr=c(0,1,0,1))
# still some bugs here:
arc.plot(x = setNames(x$Enrichment, x$Term), up.list = x$Up, down.list = x$Down,
         x.mar = c(.9,0.3), x.scale = 1.7, y.mar = c(0, 0), main = "", t = 0)
# x.mar = c(1,0.3), scale = 2
text(0, 1.1, "C", cex=2)
dev.off()

# for testing:
arc.plot(x = setNames(x$Enrichment, x$Term), up.list = x$Up, down.list = x$Down,
         x.mar = c(.5,0), x.scale = 1, y.mar = c(0, 0), main = "", t = 0)

# x <- golubGO$golub
# arc.plot( x = setNames(x$LogEnrich, x$Term), up.list = x$Upregulated, down.list = x$Downregulated, 
#           y.mar = c(0, 0.1), x.mar = c(.95, 0.2), x.scale = 2.6,
#           main = "", x.ticks = 5, x.bound = 2 )
# text(0, 1.1, "C", cex=2)

########################
### paper figure 2
########################

library(CellPlot)
library(miscset)
data(leukemiasGO)

pdf("~/tmp/fig2.pdf", height = 5)
y <- lapply(leukemiasGO, function (x) {
  x$Enrichment <- x$Significant / x$Expected
  x$Upregulated <- sapply(x$log2FoldChange, function (z) sum(z>0))
  x$Downregulated <- sapply(x$log2FoldChange, function (z) sum(z<0))
  x
})
yterms <- unique(unlist(lapply(y, function(x){
  x <- subset(x, Annotated > 5 & pvalCutOff <= 0.05 & !duplicated(ENSEMBL))
  head(sort(x, TRUE, "Enrichment"), 20)$GO.ID
})))
par(mar = c(0,.5,2.5,8))
go.histogram(y, alpha.term = "pvalCutOff", min.genes = 0, max.genes = 1e10, go.selection = yterms,
             main = "GO enrichment\nin leukemia differential gene expression\ncompared to control samples", 
             axis.cex = 1, lab.cex = 1.5)
#reorder = F,
dev.off()


########################
### TESTING
########################

############### gohist

library(CellPlot)
data(leukemiasGO)
x <- leukemiasGO
x <- lapply(x, head)
gohist(x)
gohist(x, name.label = "Term", alpha = 1e-5)
gohist(x, alpha = 1e-5, size.label = .4)

### EOF
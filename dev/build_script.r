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
### paper figure 1
########################

library(CellPlot)
data(leukemiasGO)

x <- subset(leukemiasGO$CLL, pvalCutOff <= 0.05 & Significant > 20)
x <- x[order(-x$LogEnrich),]
x$up <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i>0] })
x$dwn <- lapply(Map(setNames, x$log2FoldChange, x$GenesSignificant), function (i) { i[i<0] })

pdf("~/tmp/fig1.pdf", height = 8, width = 7)
layout(matrix(1:3,nrow=3))
cell.plot(x = setNames(x$LogEnrich, x$Term), 
          cells = x$log2FoldChange, 
          main ="GO enrichment (NoT vs CLL) and DEG directionality", 
          x.mar = c(.37, 0.01), 
          key.n = 7, 
          y.mar = c(.3, 0.3), 
          cex = 1.6, 
          cell.outer = 3, 
          bar.scale = .45, 
          space = .2)
text(0, 1.1, "A", cex=2)
sym.plot(x = setNames(x$LogEnrich, x$Term), 
         cells = x$log2FoldChange, 
         x.annotated = x$Annotated, 
         main = "",
         x.mar = c(.37, 0), y.mar = c(.45,-.5),
         key.n = 7, 
         cex = 1.6, 
         axis.cex = .8, 
         group.cex = .7)
text(0, 1.1, "B", cex=2)
arc.plot(x = setNames(x$LogEnrich, x$Term), 
         up.list = x$up, main = "",
         down.list = x$dwn, x.scale = 2.4,
         y.mar = c(0),
         x.mar = c(.9,0),
         x.bound = 2.5)
text(0, 1.1, "C", cex=2)
dev.off()

########################
### paper figure 2
########################

library(CellPlot)
data(leukemiasGO)

y <- lapply(leukemiasGO, function (x) {
  x$Upregulated <- sapply(x$log2FoldChange, function (z) sum(z>0))
  x$Downregulated <- sapply(x$log2FoldChange, function (z) sum(z<0))
  x
})
yterms <- unique(unlist(lapply(y, function(x){
  x <- subset(x, pvalCutOff <= 0.05)
  x <- x[order(x$LogEnrich),]
  head(x, 9)$GO.ID
})))

pdf("~/tmp/fig2.pdf", height = 8, width = 7)
par(mar = c(0,.5,2.5,8))
go.histogram(y, go.alpha.term = "pvalCutOff", gene.alpha.term = "padj", 
             min.genes = 5, max.genes = 1e10, go.selection = yterms, show.ttest = T,
             main = "GO enrichment\nin leukemia differential gene expression\ncompared to control samples", 
             axis.cex = 1, lab.cex = 1.5, main.cex = 1.5)
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
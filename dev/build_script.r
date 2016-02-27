### build script for the CellPlot package

# set working directory to package source folder root

########################
### package
########################

library(devtools)
#library(rmarkdown)
#render('vignettes/CellPlotManual.Rmd')
document()
install(build_vignettes = F)
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

#source("http://bioconductor.org/biocLite.R")
#biocLite("hu6800.db")
#biocLite("multtest")
#biocLite("topGO")

library(plyr)
library(multtest)
library(annotate)
library(hu6800.db)
library(topGO)

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
}, mc.cores = 10))
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

GO <- lapply(list(golub=GO), function (x) {
  t <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  s <- getSigGroups(x, t)
  r <- GenTable(x, pvalCutOff = s, topNodes = 20) #length(x@graph@nodes)
  r$LogEnrich <- r$Significant / r$Expected
  return(r)
})

golubGO <- Map(mergeGOdeg, GO, list(DEG), list(M), map.gene = "PROBEID", deg.p = "padj",deg.lfc="log2fc")
golubGO <- lapply(golubGO, subset, !is.na(PROBEID))

save(golubGO, file = "data/golubGO.rdata")

########################
### dataset leukemiasGO
########################

#source("http://bioconductor.org/biocLite.R")
#biocLite("hu6800.db")
#biocLite("leukemiasEset")
#biocLite("topGO")

library(leukemiasEset)
library(DESeq2)
library(annotate)
library(hu6800.db)
library(plyr)
library(topGO)
library(BiocParallel)
library(parallel)
library(miscset)
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
# DEG <- estimateSizeFactors(DEG)
# DEG <- estimateDispersionsGeneEst(DEG)
# dispersions(DEG) <- mcols(DEG)$dispGeneEst
# DEG <- nbinomLRT(DEG, reduced = ~ 1)
DEG <- DESeq(DEG, fitType = "mean", parallel = T, BPPARAM = MulticoreParam(20))

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
  r <- GenTable(x, pvalCutOff = s, topNodes = length(x@graph@nodes))
  r$LogEnrich <- r$Significant / r$Expected
  return(r)
}, mc.cores = length(GO))

leukemiasGO <- Map(mergeGOdeg, GOsig, DEG, list(M), map.gene = "ENSEMBL")
leukemiasGO <- lapply(leukemiasGO, subset, !is.na(ENSEMBL))
leukemiasGO <- lapply(leukemiasGO, sort, TRUE, by = "LogEnrich")

save(leukemiasGO, file = "data/leukemiasGO.rdata")

########################
### paper figure 1
########################

library(CellPlot)
data(golubGO)
data(leukemiasGO)


#svg("~/tmp/fig1.svg", height = 10)
pdf("~/tmp/fig1.pdf", height = 10)
x <- head(subset(leukemiasGO$CLL, Annotated > 5), 8)
layout(matrix(1:3,nrow=3))
par(mar=c(3,0,4,0))
par(usr=c(0,1,0,1))
cell.plot( x = setNames(x$LogEnrich, x$Term), cells = x$log2FoldChange,  x.mar = c(0.3,0), y.mar = c(0.1,0.3),
           main = "GO enrichment in NoL vs CLL differential gene expression")
text(0, 1.1, "A", cex=2)
par(usr=c(0,1,0,1))
sym.plot( x = setNames(x$LogEnrich, x$Term), cells = x$log2FoldChange, x.annotated = x$Annotated,
          x.mar = c(0.3,0), y.mar = c(0.3, 0.1), cex = 1.6, ticksize = 5, key.lab = "Enrichment")
text(0, 1.1, "B", cex=2)
par(mar=c(3,0,4,0))
par(usr=c(0,1,0,1))
arc.plot(x = setNames(x$LogEnrich, x$Term), up.list = x$Upregulated, down.list = x$Downregulated,
         x.mar = c(.95,.2), x.scale = 2.6, x.bound = 2, y.mar = c(0, 0.1), main = "")
text(0, 1.1, "C", cex=2)
dev.off()

# x <- golubGO$golub
# arc.plot( x = setNames(x$LogEnrich, x$Term), up.list = x$Upregulated, down.list = x$Downregulated, 
#           y.mar = c(0, 0.1), x.mar = c(.95, 0.2), x.scale = 2.6,
#           main = "", x.ticks = 5, x.bound = 2 )
# text(0, 1.1, "C", cex=2)

########################
### paper figure 2
########################

pdf("~/tmp/fig2.pdf", height = 5)
xg <- head(subset(leukemiasGO$CLL, Annotated > 5), 8)
xg <- xg$GO.ID
x <- leukemiasGO
x <- lapply(leukemiasGO, function(i) subset(i, GO.ID %in% x))
x <- lapply(x, function(y) {y$Upregulated <- sapply(y$Upregulated, length);y$Downregulated <- sapply(y$Downregulated, length);y})
par(mar = c(0,.5,1.5,6))
go.histogram(x, alpha.term = "pvalCutOff", min.sig = 0, main = "", min.genes = 0, reorder = F,axis.cex = 1,lab.cex = 2,
             go.selection = xg, show.ttest = F)
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
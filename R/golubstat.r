#' @title Data with DEG and GO enrichment analysis for \code{golub} dataset
#' 
#' @usage
#' data(golubstat)
#' 
#' @details Microarray gene expression data from leukemia study of Golub et al. (1999).
#' Processed as described in Dudoit et al. (2002). Adapted from package 
#' \code{multtest}.
#' Differential gene expression was performed using a Student t-test to compare
#' the two groups. GO annotation was done with the \code{annotate} and 
#' \code{hu6800.db} packages from Bioconductor. The \code{topGO} package was 
#' utilized to perform a GO enrichment test based on Fisher's exact test. 
#' Data of the DEG and GO enrichment was merged with the provided \link{go.enrich}
#' function. The code to create the data object is shown in the examples section.
#' 
#' @format 
#' A \code{data.frame} that contains the gene ontology enrichment and differential
#' gene expression tests data.\cr\cr
#' Description of the columns:\cr
#' \code{id} - character vector with GO ids of the GO terms\cr
#' \code{term} - character vector with GO terms\cr
#' \code{annotated} - integer vector with number of GO annotated genes per term\cr
#' \code{significant} - integer vector with number of significant differential expressed genes per term\cr
#' \code{expected} - numeric vector with expected number of differential expressed genes per term\cr
#' \code{p} - numeric vector with p values from the GO enrichment test\cr
#' \code{padj} - numeric vector with adjusted p values (FDR) from the GO enrichment test\cr
#' \code{loge} - numeric vector with the logarithm of the GO term enrichment (significant/expected)\cr
#' \code{genes} - list of names of genes differentially expressed per GO term\cr
#' \code{deg.log2fc} - list of logarithm of fold changes of all genes per GO term\cr
#' \code{deg.p.adj} - list of adjusted p values of all genes per GO term\cr
#' 
#' @references
#' 
#' Golub et al. (1999), Molecular Classification of Cancer: Class Discovery and
#' Class Prediction by Gene Expression Monitoring, 
#' Science, Volume 286 no. 5439 pp. 531-537, doi:10.1126/science.286.5439.531
#' 
#' Dudoit, Fridlyand and Speed (2002), Comparison of Discrimination Methods for
#' the Classification of Tumors Using Gene Expression Data, 
#' Journal of the American Statistical Association, Volume 97, Issue 457,
#' doi:10.1198/016214502753479248
#' 
#' @author 
#' Sven E. Templer [aut]\cr
#' Authors of multtest package [ctb]
#' 
#' @examples
#' 
#' \dontrun{
#' 
#'  ### Differential gene expression
#'  ### and GO enrichment analysis
#'  ### of the golub dataset:
#'  
#'  ## install dependencies
#'  source("http://bioconductor.org/biocLite.R")
#'  biocLite("hu6800.db")
#'  biocLite("multtest")
#'  biocLite("topGO")
#'  
#'  ## load dependencies
#'  library(cellplot)
#'  library(plyr)
#'  library(multtest)
#'  library(annotate)
#'  library(hu6800.db)
#'  library(topGO)
#'  
#'  ## data import
#'  data(golub)
#'  
#'  ## differential gene expression testing
#'  golub.test <- function (e, i, n) {
#'    i <- as.logical(i)
#'    a <- split(e[,i],seq(nrow(e)))
#'    b <- split(e[,!i],seq(nrow(e)))
#'    fun.map <- function (a, b) {
#'      t <- t.test(a, b)
#'      t$fc <- t$estimate[1]/t$estimate[2]
#'      suppressWarnings(t$fc <- log2(t$fc))
#'      data.frame(
#'        mean.aml = t$estimate[1],
#'        mean.all = t$estimate[2],
#'        log2fc  = t$fc,
#'        p = t$p.value,
#'        row.names = NULL)
#'    }
#'    s <- Map(fun.map, a, b)
#'    s <- do.call(rbind, s)
#'    s$p.adj <- p.adjust(s$p)
#'    s <- data.frame(n, s, stringsAsFactors = FALSE)
#'    s <- subset(s, !is.na(log2fc))
#'    return(s)
#'  }
#'  
#'  degdata <- golub.test(golub, golub.cl, data.frame(
#'    gene = golub.gnames[,3], index = as.integer(golub.gnames[,1]),
#'    stringsAsFactors = FALSE))
#'  
#'  ## gene ontology annotation
#'  goterms <- select(
#'    hu6800.db, degdata$gene,
#'    c("PROBEID","ALIAS","GO","ENSEMBL","ENTREZID"))
#'  goterms <- dlply(
#'    subset(goterms, ONTOLOGY == "BP"), "PROBEID",
#'    function (x) unique(x$GO))
#'  
#'  ## topGO object
#'  godata <- new(
#'    "topGOdata", ontology = "BP", description = 'golub',
#'    allGenes = setNames(degdata$p.adj, degdata$gene),
#'    geneSelectionFun = function (allScore) { allScore <= 0.05 },
#'    annotationFun = annFUN.gene2GO, gene2GO = goterms)
#'  
#'  ## gene ontology enrichment
#'  #golubstat <- go.enrich(godata, degdata, p.max = 1)
#'  golubstat <- go.enrich(godata, degdata)
#'
#  
# # store
# save(golubstat, file = "data/golubstat.rdata")
# str(golubstat, list.len = 12)
#
#' }
#' 

"golubstat"
golubstat <- NULL

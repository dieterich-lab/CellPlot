#' @title Data with DEG and GO enrichment analysis for \code{golub} dataset
#' 
#' @details Microarray gene expression data from leukemia study of Golub et al. (1999).
#' Processed as described in Dudoit et al. (2002). Adapted from package 
#' \code{multtest}.
#' Differential gene expression was performed using a Student t-test to compare
#' the two groups. See details in \link{golub.degtest}.
#' GO annotation has been performed with the \code{annotate} and \code{hu6800.db}
#' packages from Bioconductor. The \code{topGO} package was used to perform
#' a GO enrichment test based on Fisher's exact test. Data of the DEG and GO
#' enrichment was merged with the provided \link{topgo2cellplot} function.
#' The code to create the data object is shown in the examples section.
#' 
#' @format A list containing 2 slots:
#' \describe{
#'   \item{\code{stats}}{A data.frame with the gene expression statistics.}
#'   \item{\code{go}}{A list with the merged DEG and GO enrichment results.}
#' }
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
#' \dontrun{
#' 
#' ### Differential gene expression
#' ### and GO enrichment analysis:
#' 
#' ## install dependencies
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("hu6800.db")
#' biocLite("multtest")
#' biocLite("topGO")
#' 
#' ## load dependencies
#' library(plyr)
#' library(multtest)
#' library(annotate)
#' library(hu6800.db)
#' library(topGO)
#' 
#' ## data import and gene golub test
#' data(golub)
#' golub.deg <- list()
#' golub.deg$stats <- golub.degtest(
#'   golub, as.logical(golub.cl), 
#'     data.frame(gene = golub.gnames[,3], 
#'                gene.index = as.integer(golub.gnames[,1]),
#'                stringsAsFactors = FALSE))
#'                
#' ## GO annotation
#' golub.deg$go <- select(hu6800.db, golub.deg$stats$gene, 
#'                        c("PROBEID","ALIAS","GO","ENSEMBL","ENTREZID"))
#' golub.deg$go <- dlply(subset(golub.deg$go, ONTOLOGY == "BP"), "PROBEID", 
#'                       function (a) unique(a$GO))
#'                       
#' ## build topGO object with GO topology
#' golub.deg$go <- new(
#'   "topGOdata", ontology="BP", description='golub.deg',
#'   allGenes=setNames(golub.deg$stats$p.adj, golub.deg$stats$gene),
#'   geneSelectionFun=function(allScore){allScore <= 0.05},
#'   annotationFun=annFUN.gene2GO, gene2GO=golub.deg$go)
#'   
#' ## GO enrichment test and data merge with DEG results
#' golub.deg$go <- topgo2cellplot(golub.deg$go, golub.deg$stats$log2fc, 
#'                                golub.deg$stats$gene)
#'
#' ## store data
#' # save(golub.deg, file = "data/golub.deg.rdata")
#' # str(golub.deg$go, list.len = 8)
#' }
#' 

"golub.deg"

golub.deg <- NULL

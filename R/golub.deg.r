#' @title DEG and GO enrichment of \code{golub} data
#' 
#' @details Microarray gene expression data from leukemia study of Golub et al. (1999).
#' Processed as described in Dudoit et al. (2002). Adapted from package 
#' \code{multtest}.
#' Differential gene expression analysis and GO annotation and enrichment was 
#' performed as shown in the example section.
#' 
#' @format A list containing 2 slots:
#' \describe{
#'   \item{stats}{A data.frame with the gene expression statistics.}
#'   \item{go}{A list with the GO enrichment statistics.}
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
#' @author Sven E. Templer
#' 
#' @examples
#' \dontrun{
#' 
#' ### Differential gene expression
#' ### and GO enrichment analysis:
#' 
#' ## install dependencies
#' install.packages("multtest")
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("hu6800.db")
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
#' golub.deg$stats <- golub.deg.test(
#'   golub, as.logical(golub.cl), 
#'     data.frame(gene = golub.gnames[,3], 
#'                gene.index = as.integer(golub.gnames[,1]),
#'                stringsAsFactors = FALSE))
#'                
#' ## annotation and GO enrichment test
#' golub.deg$go <- select(hu6800.db, golub.deg$stats$gene, 
#'                        c("PROBEID","ALIAS","GO","ENSEMBL","ENTREZID"))
#' golub.deg$go <- dlply(subset(golub.deg$go, ONTOLOGY == "BP"), "PROBEID", 
#'                       function (a) unique(a$GO))
#' golub.deg$go <- new(
#'   "topGOdata", ontology="BP", description='golub.deg',
#'   allGenes=setNames(golub.deg$stats$p.adj, golub.deg$stats$gene),
#'   geneSelectionFun=function(allScore){allScore <= 0.05},
#'   annotationFun=annFUN.gene2GO, gene2GO=golub.deg$go)
#' golub.deg$go <- topgo2cellplot(golub.deg$go, golub.deg$stats$log2fc, 
#'                                golub.deg$stats$gene)
#' 
#' ### Visualization with cellplot:
#' 
#' cell.plot(golub.deg$go$go.loge, golub.deg$go$deg.logfc)
#' 
#' }
#' 

"golub.deg"

#save(golub.deg, file = "data/golub.deg.rdata")
#str(golub.deg$go, list.len = 5)

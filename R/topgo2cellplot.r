#' @title GO enrichment test and conversion of topGO data to list for cellplot
#' 
#' @description The function performs the GO enrichment tests (a Fisher test).
#' It extracts data from \code{topGO} objects and combines it with
#' differential gene expression fold changes to create a data list
#' suitable for \code{cellplot} input.
#' 
#' @param go A \code{topGO} object.
#' @param deg.logfc A numeric vector with the log fold changes from the
#' differential gene expression analysis.
#' @param deg.genes A character vector with the gene names corresponding
#' to \code{deg.logfc}.
#' @param n An integer value referring to the number of top hits to extract.
#' @param ids A logical value to attach GO IDs to the GO terms.
#' @param replace.infinite A logical value to replace infinite values
#' in \code{deg.logfc} with absolute maximums plus/minus one, depending on
#' the sign of the maximums.
#' @param p.adj.method A character string naming the p value adjustment
#' method. Use \code{"none"} to pass through. See \link{p.adjust}.
#' 
#' @author Sven E. Templer

#' @export
topgo2cellplot <- function (
  go, deg.logfc = NULL, deg.genes = NULL, n = 30, ids = FALSE, 
  replace.infinite = FALSE, p.adj.method = 'fdr')
{
  # go: topGOdata
  # cuf: cuffdiff # x <- read.cuffdiff()
  # n: number of top nodes
  
  # test GO enrichment
  go.test <- getSigGroups(go, new("classicCount", testStatistic = GOFisherTest, name = "Fisher test"))
  go.table <- GenTable(go, p=go.test, orderBy="p", ranksOf="p", topNodes=n)
  go.table$p <- as.numeric(go.table$p)
  gn <- length(go@graph@nodes)
  #gn <- length(topGO:::genes(go))
  ret.logp <- -log2(p.adjust(go.table$p, method = p.adj.method, n=gn))
  ret.loge <- log(go.table$Significant / go.table$Expected)
  
  # convert data
  go.signames <- sigGenes(go)
  ret.genes <- genesInTerm(go, go.table$GO.ID)
  ret.genes <- lapply(ret.genes, function (i) {intersect(i,go.signames)})
  ret.term <- go.table$Term
  ret.ids <- names(ret.genes)
  
  # check if deg.logfc
  if (!is.null(deg.logfc)) {
    c.logfc <- setNames(deg.logfc, deg.genes)
    ret.logfc <- lapply(ret.genes, function (i) { setNames(c.logfc[i], NULL) })
    if (replace.infinite) {
      warning("replace.infinite will be deprecated")
      ret.logfc <- lapply(ret.logfc, function(i) {
        inf <- is.infinite(i)
        i[inf & i>0] <- max(i[!inf],0) + 1
        i[inf & i<0] <- min(i[!inf],0) - 1
        i
      })
    } #else { warning("replace.infinite=TRUE necessary for cell.plot!") }
  } else {
    ret.logfc <- lapply(ret.genes, function (i) rep(.0, length(i)))
  }
  
  ret.n <- setNames(sapply(ret.logfc, length), NULL)
  o <- order(ret.loge, decreasing = T)
  ret <- list(
    deg.logp  = ret.logp[o],
    deg.logfc = ret.logfc[o],
    go.loge   = ret.loge[o],
    go.ids    = ret.ids[o],
    go.terms  = ret.term[o],
    n         = ret.n[o],
    genes     = ret.genes[o])
  id.names <- if (ids) paste(ret$go.terms, ret$go.ids) else ret$go.terms
  ret$go.loge <- setNames(ret$go.loge, id.names)
  ret$deg.logp <- setNames(ret$deg.logp, id.names)
  attr(ret, "gotab") <- go.table
  return(ret)
}
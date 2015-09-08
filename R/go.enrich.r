#' @title GO enrichment test and conversion of topGO data to list for cellplot
#' 
#' @description The function performs the GO enrichment tests and merges
#' the gene expression statistics into the output. A \code{topGO} object
#' needs to be provided as input to run the enrichment test on, as well
#' as a corresponding data.frame containing the gene expression tests.
#' 
#' @param go A \code{topGO} object.
#' @param deg A data frame with the differential gene expression statistics
#' and gene identifiers.
#' @param col.stat A character vector with the column name(s) from \code{deg}
#' that contain statistics, e.g. log fold change or (adjusted) p values.
#' @param col.genes A character vector with the column name containing the
#' gene identifiers.
#' @param n.top An integer value referring to the number of top hits to extract.
#' @param test A \code{topGO} test function, e.g. \link[topGO]{GOFisherTest}.
#' Used as \code{testStatistic} argument value in the \link[topGO]{classicCount-class}
#' creator function.
#' @param name A character string for the test.
#' Used as \code{name} argument value in the \link[topGO]{classicCount-class}
#' creator function.
#' @param p.max A numeric value for the significance threshold (p value)
#' of the GO enrichment test. Set to \code{1} to extract all results from
#' the test.
#' @param p.adj.method A character string naming the p value adjustment
#' method. Use \code{"none"} to pass through. See \link{p.adjust}.
#' @param order.by A character vector, either \code{"e"} or \code{"p"} to sort
#' the output by enrichment or p values of the GO enrichment analysis.
#' 
#' @return 
#' The result is a \code{data.frame}. It contains the output of \code{GenTable}
#' from the \code{topGO} package merged with gene names and expression statistics
#' from a differential gene expression data frame.\cr
#' Columns are described at the \code{\link{golubstat}} data help page. Find an
#' example workflow in the examples section, too.\cr
#' The output is suitable for the \link{cell.plot} or \link{sym.plot} functions.\cr
#' 
#' @author 
#' Sven E. Templer [aut]
#' 
#' @name go.enrich

NULL

GOFisherTest <- topGO::GOFisherTest
#GOFisherTest <- NULL

#' @rdname go.enrich
#' @export 
go.enrich <- function (
  go, deg, col.stats = c("log2fc","p.adj"), col.genes = "gene", 
  n.top = NULL, test = GOFisherTest, name = "Fisher test",
  p.max = 0.05, p.adj.method = 'fdr', order.by = c("enrichment","pvalue"))
{
  
  # checks
  if (!requireNamespace("topGO")) 
    stop("Install package topGO")
  order.by <- match.arg(order.by)
  if (!all(col.stats %in% names(deg))) 
    stop("Select valid stats columns from deg")
  if (length(col.genes) != 1 || !col.genes %in% names(deg)) 
    stop("Select valid genes columns from deg")
  
  # GO enrichment
  go.n <- length(go@graph@nodes)
  go.test <- new("classicCount", testStatistic = test, name = name)
  go.test <- topGO::getSigGroups(go, go.test)
  res <- topGO::GenTable(go, p = go.test, topNodes = go.n)
  names(res)[1:5] <- c("id", "term", "annotated", "significant", "expected")
  res$p <- as.numeric(res$p)
  res$padj <- p.adjust(res$p, method = p.adj.method, n = go.n)
  res$loge <- log(res$significant / res$expected)
  
  # add deg
  go.signames <- topGO::sigGenes(go)
  go.genes <- topGO::genesInTerm(go, res$id)
  go.genes <- lapply(go.genes, function (i) { intersect(i, go.signames) })
  deg.stats <- deg[, col.stats, drop = F]
  deg.stats <- lapply(deg.stats, setNames, nm = deg[[col.genes]])
  deg.genes <- lapply(go.genes, function (g) { names(deg.stats[[1]][g]) })
  names(deg.genes) <- NULL
  deg.stats <- lapply(deg.stats, function (s) {
    lapply(go.genes, function (g) { setNames(s[g], NULL) })
  })
  res$genes <- deg.genes
  names(deg.stats) <- paste0("deg.", names(deg.stats))
  for (name in names(deg.stats))
    res[[name]] <- setNames(deg.stats[[name]], NULL)
  
  # order
  res <- subset(res, p <= p.max)
  o <- order(switch(order.by, enrichment = abs(res$loge), pvalue = -log(res$p)), decreasing = T)
  res <- res[o,]
  if (is.null(n.top) || n.top > nrow(res))
    n.top <- nrow(res)
  res <- res[seq(n.top),]
  rownames(res) <- NULL
  return(res)
  
}

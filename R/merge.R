#' @title Merge topGO results with gene expression tables
#' 
#' @description 
#' Merge a topGO result from the \code{GenTable} output with a gene expression table
#' to have a list of gene names and statistics for each GO term.
#' 
#' @param go The topGO \code{GenTable()} data.frame.
#' @param deg The gene expression data.frame. See details.
#' @param map The mapping table for GO term to gene name conversion. See details.
#' @param go.go The column name for the GO identifiers.
#' @param deg.gene The column name for the gene identifiers.
#' @param deg.stats The column name(s) for the gene expression statistics
#' (e.g. p-values or log2 fold changes).
#' @param map.go The column name for the GO identifiers.
#' @param map.gene The column name for the gene identifiers.
#' 
#' @details 
#' tba
#' 
#' @author 
#' Sven E. Templer [aut]

#' @export
mergeGOdeg <- function (
  go, deg, map, 
  go.go = "GO.ID",
  deg.gene = "gene", deg.p = "pvalue", deg.lfc = "log2FoldChange",
  map.go = "GO", map.gene = "gene",
  as.intersect = FALSE)
{
  
  if (!go.go %in% colnames(go)) stop("missing go.go")
  if (!deg.gene %in% colnames(deg)) stop("missing deg.gene")
  if (!deg.p %in% colnames(deg)) stop("missing deg.p")
  if (!is.null(deg.lfc) && !deg.lfc %in% colnames(deg)) stop("missing deg.lfc")
  if (!map.go %in% colnames(map)) stop("missing map.go")
  if (!map.gene %in% colnames(map)) stop("missing map.gene")
  
  cat("* counts\n",
      " go.go               ", length(unique(go[[go.go]])), "\n",
      " map.go              ", length(unique(map[[map.go]])), "\n",
      " deg.gene            ", length(unique(deg[[deg.gene]])), "\n",
      " map.gene            ", length(unique(map[[map.gene]])), "\n",
      " go.go & map.go      ", length(unique(intersect(go[[go.go]],map[[map.go]]))), "\n",
      " deg.gene & map.gene ", length(unique(intersect(deg[[deg.gene]],map[[map.gene]]))), "\n"
      )
  if (as.intersect) {
    isect.go <- intersect(go[[go.go]], map[[map.go]])
    isect.gene <- intersect(deg[[deg.gene]], map[[map.gene]])
    map.isect <- 
      map[[map.go]] %in% isect.go &
      map[[map.gene]] %in% isect.gene
    map <- map[map.isect,]
    cat("  intersect go        ", length(unique(map[[map.go]])), "\n",
        " intersect deg       ", length(unique(map[[map.gene]])), "\n")
    if (!nrow(map))
      stop("mapping resulted in no intersect")
    go <- go[go[[go.go]] %in% map[[map.go]],]
    deg <- deg[deg[[deg.gene]] %in% map[[map.gene]],]
  }
  
  x <- merge(map[,c(map.go, map.gene)], deg[,c(deg.gene, deg.p, deg.lfc)],
             by.x = map.gene, by.y = deg.gene) # all = TRUE?
  x <- merge(go, x, by.x = go.go, by.y = map.go, all.x = TRUE, all.y = FALSE)
  gi <- match(go.go, colnames(x))
  x <- aggregate(x[,-gi], by=list(x[[gi]]), function (i) if (length(unique(i))>1) i else unique(i))
  colnames(x)[1] <- go.go
  if (length(x$pvalCutOff)) x$pvalCutOff <- as.numeric(x$pvalCutOff)
  if (!is.null(deg.lfc)) {
    xfc <- Map(setNames, x[[deg.lfc]], x[[map.gene]])
    x$Upregulated <- lapply(xfc, function (i) { i[i>0] })
    x$Downregulated <- lapply(xfc, function (i) { i[i<0] })
  }
  return(x)
}


#' @title Perform DEG and GO enrichment analysis (of \code{golub} data)
#' 
#' @description  Method used to analyse differential gene expression of the
#' \code{golub} data to create the \code{golub.deg} data.
#' 
#' @param e Matrix with expression values as in \code{golub} data.
#' Genes per rows, samples per column.
#' @param i Logical vector as column/sample index for \code{e}, for
#' example the leukemia classes from \code{golub} data.
#' @param a.frame A \code{data.frame} which can contains annotation information
#' like gene names and which is included into the output.
#' 
#' @return Returns a \code{data.frame}. See \link{golub.deg} for more
#' information.
#' 
#' @author 
#' Sven E. Templer [aut]

#' @export
golub.degtest <- function (e, i, a.frame = NULL) {
  s <- Map(function (a,b) {
    t <- t.test(a,b)
    t$fc <- t$estimate[1]/t$estimate[2]
    suppressWarnings(t$fc <- log2(t$fc))
    data.frame(
      mean.aml = t$estimate[1],
      mean.all = t$estimate[2],
      log2fc  = t$fc,
      p = t$p.value,
      row.names = NULL)
  }, a= split(e[,i],seq(nrow(e))), b = split(e[,!i],seq(nrow(e))))
  s <- do.call(rbind, s)
  s$p.adj <- p.adjust(s$p)
  if (!is.null(a.frame)) 
    s <- data.frame(a.frame, s, stringsAsFactors = FALSE)
  s <- subset(s, !is.na(log2fc))
  return(s)
}

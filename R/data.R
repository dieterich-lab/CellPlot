#' @name leukemiasGO
#' @aliases golubGO
#' @title Gene ontology term enrichment of differential gene expression of leukemia samples
#' 
#' @description 
#' \code{golubGO} is microarray gene expression data from leukemia study of 
#' Golub et al. (1999). Processed as described in Dudoit et al. (2002). Adapted
#' from package \code{multtest} (see \link[multtest]{golub}). Differential gene
#' expression was performed using a t-test to compare the two groups ALL and
#' AML.\cr
#' \cr
#' \code{leukemiasGO} is microarray gene expression data from the studies of 
#' Kohlmann et al. (2008) and Haferlach et al. (2010). Processed as described
#' in Risueno et al. (2010). Adapted from package \code{leukemiasEset} (see
#' \link[leukemiasEset]{leukemiasEset}). Differential gene expression was 
#' performed using the \code{DESeq2} package comparing the
#' sample classes NoL (no leukemia) against AML, ALL, and CML.\cr
#' \cr
#' GO annotation was done with the \code{annotate} and \code{hu6800.db} 
#' packages from Bioconductor. The \code{topGO} package was utilized to perform
#' a GO enrichment test based on Fisher's exact test using the 'elimScore' 
#' method.
#' 
#' @usage
#' data(golubGO)
#' data(leukemiasGO)
#' 
#' @format 
#' A \code{data.frame} that contains the gene ontology enrichment and differential
#' gene expression tests data.\cr\cr
#' Description of the columns:\cr
#' \code{GO} - character vector with GO ids of the GO terms\cr
#' \code{Term} - character vector with GO terms\cr
#' \code{Annotated} - integer vector with number of GO annotated genes per term\cr
#' \code{Significant} - integer vector with number of significant differential expressed genes per term\cr
#' \code{Expected} - numeric vector with expected number of differential expressed genes per term\cr
#' \code{pvalCutOff} - numeric vector with p values from the GO enrichment test\cr
#' \code{padj} - numeric vector with adjusted p values (FDR) from the GO enrichment test\cr
#' \code{LogEnrich} - numeric vector with the logarithm of the GO term enrichment (significant/expected)\cr
#' \code{PROBEID/ENSEMBL} - list of names of genes differentially expressed per GO term\cr
#' \code{log2fc/log2FoldChange} - list of logarithm of fold changes of all genes per GO term\cr
#' \code{p/pvalue} - list of p values of all differential gene expression\cr
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
#' Kohlmann A, Kipps TJ, Rassenti LZ, Downing JR et al. An international 
#' standardization programme towards the application of gene expression 
#' profiling in routine leukaemia diagnostics: the Microarray Innovations in 
#' LEukemia study prephase. Br J Haematol (2008) 142(5):802-7. PMID: 18573112
#' 
#' Haferlach T, Kohlmann A, Wieczorek L, Basso G et al. Clinical utility of 
#' microarray-based gene expression profiling in the diagnosis and 
#' subclassification of leukemia: report from the International Microarray 
#' Innovations in Leukemia Study Group. J Clin Oncol (2010) 28(15):2529-37. 
#' PMID: 20406941
#' 
#' Risueno A, Fontanillo C, Dinger ME, De Las Rivas J. GATExplorer: genomic and
#' transcriptomic explorer; mapping expression probes to gene loci, transcripts,
#' exons and ncRNAs. BMC Bioinformatics (2010) 11:221. PMID: 20429936. 
#' 
#' @author 
#' Sven E. Templer [aut]\cr
#' Authors of the multtest package [ctb]\cr
#' Authors of the leukemiasEset [ctb]
#' 

"golubGO"
golubGO <- NULL

"leukemiasGO"
leukemiasGO <- NULL
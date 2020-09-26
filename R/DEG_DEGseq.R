#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @noRd


.DEG_DEGseq <- function(countTable,
                      FDR_threshold = .1){

    ## Identify DEG with edgeR
    ##
    ## Args:
    ##   countTable: a count matrix of the RNASeq data
    ##   FDR_threshold: a threshold used by DEGseq to determine DEG
    ##
    ## Return:
    ##  A data.frame of various statistics regarding differential expression.
    ##  Each row corresponds to a gene.


    m1 <- as.matrix(cbind( gene_name = rownames(countTable), countTable[1]))
    m2 <- as.matrix(cbind( gene_name = rownames(countTable), countTable[2]))
    folder = paste0('./', sample(1:1000000,1))
    DEGseq::DEGexp(geneExpMatrix1 = m1,
           geneCol1 = 1,
           expCol1 = 2,
           groupLabel1 = colnames(m1)[2],
           geneExpMatrix2 = m2,
           geneCol2 = 1,
           expCol2 = 2,
           groupLabel2 = colnames(m2)[2],
           method = 'MARS',
           thresholdKind = 3, #BH FDR
           qValue = FDR_threshold,
           outputDir = folder)
    res <- utils::read.delim(file = paste0(folder,'/output_score.txt'))
    unlink(folder, recursive = T)
    return(res)
}

#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#' @noRd
#'
.DEG_DEGseq.cb <- function(countTable1,countTable2,
                       FDR_threshold = .1){
  ## Identify DEG with edgeR
  ##
  ## Args:
  ##   countTable: a count matrix of the RNASeq data
  ##   FDR_threshold: a threshold used by DEGseq to determine DEG
  ##
  ## Return:
  ##  A data.frame of various statistics regarding differential expression.
  ##  Each row corresponds to a gene.


  m1 <- as.matrix(cbind( gene_name = rownames(countTable1), countTable1))
  m2 <- as.matrix(cbind( gene_name = rownames(countTable2), countTable2))
  folder = paste0('./', sample(1:1000000,1))
  DEGseq::DEGexp(geneExpMatrix1 = m1,
         geneCol1 = 1,
         expCol1 = 2,
         groupLabel1 = colnames(m1)[2],
         geneExpMatrix2 = m2,
         geneCol2 = 1,
         expCol2 = 2,
         groupLabel2 = colnames(m2)[2],
         method = 'MARS',
         thresholdKind = 3, #BH FDR
         qValue = FDR_threshold,
         outputDir = folder)
  res <- utils::read.delim(file = paste0(folder,'/output_score.txt'))
  unlink(folder, recursive = T)
  return(res)
}

#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @noRd


.DEG_DESeq2 <- function(countTable,conditions){


    ## Identify DEG w/ DEGSeq2
    ##
    ## Args:
    ##   countTable: a count matrix of the RNASeq data
    ##   condtions: the experiment conditions of the two samples
    ##
    ## Returns:
    ##  A data.frame of statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    ## produce colData required by DESeq2,
    ## and make the format conform with what DESeq2 requires
    colDataTmp <- data.frame(condition = conditions)
    row.names(colDataTmp) <- conditions
    ## convert count data matrix to the formated that DESeq operates on
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countTable,
                                  colData =colDataTmp,
                                  design = ~condition)
    ## conduct DESeq2
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    return(res)
}

.DEG_DESeq2.cb <- function(countTable,conditions){
    ## Identify DEG w/ DEGSeq2
    ##
    ## Args:
    ##   countTable: a count matrix of the RNASeq data
    ##   condtions: the experiment conditions of the two samples
    ##
    ## Returns:
    ##  A data.frame of statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    ## produce colData required by DESeq2,
    ## and make the format conform with what DESeq2 requires
    colDataTmp <- data.frame(condition = conditions)
    #row.names(colDataTmp) <- conditions
    ## convert count data matrix to the formated that DESeq operates on
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countTable,
                                  colData =colDataTmp,
                                  design = ~condition)
    ## conduct DESeq2
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    return(res)
}



#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @noRd


.DEG_DESeq <- function(countTable,conditions){

  # require('locfit')
  ## Identify DEG w/ DESeq
    ##
    ## Args:
    ##  countTable: a count matrix of the RNASeq data
    ##  condtions: the experiment conditions of the two samples
    ##
    ## Returns:
    ##  A data.frame of statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    designDf = data.frame(row.names = c('case_rep1','control_rep1'),
                          condition = conditions)
    cds2 = DESeq::newCountDataSet(countTable, designDf$condition)
    cds2 = DESeq::estimateSizeFactors(cds2)
    cds2_try = try(DESeq::estimateDispersions(cds2, method='blind', sharingMode = 'fit-only'),silent = T) # the only step that is different for analysis w/o replicates
    if(class(cds2_try) == 'try-error'){
        message('Parametric dispersion fit failed, tried a local fit. ')
        cds2_try = DESeq::estimateDispersions(cds2, method='blind', sharingMode = 'fit-only', fitType='local')
    }
    res = DESeq::nbinomTest(cds2_try,colnames(countTable)[1],colnames(countTable)[2])
    return(res)
}

####### Cohort Based DESeq

.DEG_DESeq.cb <- function(countTable,conditions){
  # require('locfit')

    ## Identify DEG w/ DESeq
    ##
    ## Args:
    ##  countTable: a count matrix of the RNASeq data
    ##  condtions: the experiment conditions of the two samples
    ##
    ## Returns:
    ##  A data.frame of statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    #  designDf = data.frame(row.names = c('case','control', 'case','control','case','control'),
    #                      condition = conditions)

    cds2 = DESeq::newCountDataSet(countTable, conditions)
    cds2 = DESeq::estimateSizeFactors(cds2)
    cds2_try = try(DESeq::estimateDispersions(cds2, method='per-condition', sharingMode = 'fit-only'),silent = T) # the only step that is different for analysis w/o replicates
    if(class(cds2_try) == 'try-error'){
        message('Parametric dispersion fit failed, tried a local fit. ')
        cds2_try = DESeq::estimateDispersions(cds2, method='per-condition', sharingMode = 'fit-only', fitType='local')
    }
    res = DESeq::nbinomTest(cds2_try,unique(conditions)[1],unique(conditions)[2])
    return(res)
}

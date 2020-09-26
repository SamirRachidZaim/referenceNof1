#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @noRd



.DEG_edgeR <- function(countTable,
                      conditions,
                      replicateType = c('technical replicates',
                                        'genetically identical model organisms',
                                        'human data','user specified'), disper){

  ## Identify DEG with edgeR
    ##
    ## Args:
    ##   countTable: a count matrix of the RNASeq data.
    ##               Note this matrix has only two columns.
    ##   condtions: the experiment conditions of the two samples
    ##   disper: a user provided estimate of the common value of the dispersion
    ##           parameter shared by all genes.
    ##   replicateType: types of replicates (see manual of edgeR for details)
    ##
    ## Return:
    ##  A data.frame of various statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    y <- edgeR::DGEList(countTable, group = factor(conditions))
    y <- edgeR::calcNormFactors(y)
    repType = match.arg(replicateType)
    if (repType == 'technical replicates'){
        bcv = .001
    }else if(repType == 'genetically identical model organisms'){
        bcv = .1
    }else if(repType == 'human data'){
        bcv = .4
    }else if(repType == 'user specify'){
        bcv = sqrt(disper)
    }
    et <- edgeR::exactTest(y,dispersion = bcv^2)
    return(et)
}

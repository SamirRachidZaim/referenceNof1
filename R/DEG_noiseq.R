#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @noRd


.DEG_noiseq <- function(countTable,conditions){


    ## Identify DEG w/ NOISeq
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
  mydata <- NOISeq::readData(data = countTable, factors = designDf)
  myresults <- NOISeq::noiseq(mydata, factor = "condition", k = NULL, norm = "n",
                      nss = 3, v = 0.02, replicates = "no")
  res = myresults@results[[1]]
  return(res)
}


.DEG_noiseq.cb <- function(countTable,conditions){
    ## Identify DEG w/ NOISeq
    ##
    ## Args:
    ##  countTable: a count matrix of the RNASeq data
    ##  condtions: the experiment conditions of the two samples
    ##
    ## Returns:
    ##  A data.frame of statistics regarding differential expression.
    ##  Each row corresponds to a gene.

    designDf = data.frame(row.names = colnames(countTable),
                          condition = conditions)
    mydata <- NOISeq::readData(data = countTable, factors = designDf)
    myresults = NOISeq::noiseqbio(mydata, k = 0.5, norm = "rpkm", nclust = 15, plot = FALSE, factor='condition')

    res = myresults@results[[1]]
    return(res)
}

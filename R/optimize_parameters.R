#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage optimize_parameters(expressionCutoffs = c(30,50),
#'                            FCRegions = list(c(0,0.3)),
#'                            gs.mat=gs.mat,
#'                            FDR.cutoff = 0.1)
#'
#' @param expressionCutoffs a vector of positive integer > 1 indicating the minimum count value of RNA-seq
#' @param FCRegions a tlist of tuples indicating the fold change region to search
#' @param gs.mat the matrix used to develop the reference standard
#' @param fdr.method fdr.threshold for determining which set of genes are DEGs in the reference standard
#'

optimize_parameters <- function(expressionCutoffs, FCRegions,gs.mat, FDR.cutoff=0.1){

  params = expand.grid(expressionCutoffs, FCRegions)
  colnames(params) <- c("Expression Cutoff", 'Fold Change Window')

  search_space <- lapply(1:nrow(params),
                     function(x) applyFilter(expressionCutoff =  params[x,1],
                                                            fold_change_window=unlist(params[x,2]),
                                                            gs.mat=gs.mat,
                                                            FDR.cutoff =FDR.cutoff))


 medianOfMedians <- sapply(search_space, function(x) median(x$Medians) )
 idx.max <- which.max(medianOfMedians)

 optimal.param <- params[idx.max,]
 optimal.score <- medianOfMedians[idx.max]
 optimal.param$MedianAgreement <- optimal.score


 cat('\n\nThe best Jaccard index concordances were attained using \n\n')
 cat(paste('an expression level cutoff of ', optimal.param[1]))
 cat('\n\n')
 cat(paste('and a fold-change window of ', optimal.param[2]))
 cat(paste('attaining a median JI concordance across all methods of', round(optimal.score,2)))
 cat('\n\n')
 return(optimal.param)
}
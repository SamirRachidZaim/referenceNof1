#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage optimize_parameters(gs.mat=utils::data("mcf7"),
#'                            expressionCutoffs = c(30,50),
#'                            FCRegions = list(c(0,0.3)),
#'                            FDR.cutoff = 0.1)
#'
#' @param gs.mat the matrix used to develop the reference standard
#' @param expressionCutoffs a vector of positive integer > 1 indicating the minimum count value of RNA-seq
#' @param FCCutoff a list of upper bounds for minimal fold change thresholds
#' @param FDR.cutoff fdr.threshold for determining which set of genes are DEGs in the reference standard
#' @param targetJI target jaccard index 
#'
#' @export

optimize_parameters <- function(gs.mat=utils::data("mcf7"), expressionCutoffs= c(30,50), FCCutoff=c(1,1.3,1.5), FDR.cutoff=0.1, targetJI=0.6){

  FCRegions = lapply(FCCutoff, function(x) c(x, Inf))
  
  requireNamespace('locfit')
  params = expand.grid(expressionCutoffs, FCRegions)
  colnames(params) <- c("Expression Cutoff", 'Fold Change Window')
  
  i=1
  currentJI=0
  while(currentJI < targetJI ){
    
    curr_search = applyFilter(gs.mat=gs.mat,
               expressionCutoff =  params[i,1],
               fold_change_window=unlist(params[i,2]),
               FDR.cutoff =FDR.cutoff)    
    
    currentJI = median(curr_search$Medians)
    i=i+1
  }

  optimal.idx = i-1

  curr_search

 optimal.param <- params[optimal.idx,]
 optimal.score <- round(currentJI,2)
 optimal.param$MedianAgreement <- optimal.score

 cat('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nThe best Jaccard index concordances were attained using')
 cat(paste('an expression level cutoff of ', optimal.param[1]))
 cat(paste('and a fold-change window of ', optimal.param[2]))
 cat(paste(' attaining a median JI concordance across all methods of', round(optimal.score,2)))
 cat('\n\n')
 return(optimal.param)
}

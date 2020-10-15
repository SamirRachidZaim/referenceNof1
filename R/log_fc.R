#' random forest feature selection based on binomial exact test
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage log_fc(pair.mat)
#'
#' @param pair.mat a two-column, stimulus-response paired sample matrix
#'
#' @return a data.frame with 1 columns containing the log fold change


log_fc <- function(pair.mat){
  log2(pair.mat[,2]+1)/log2(pair.mat[,1]+1)
}

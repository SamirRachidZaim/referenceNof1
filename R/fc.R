#' random forest feature selection based on binomial exact test
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage fc(pair.mat)
#'
#' @param pair.mat a two-column, stimulus-response paired sample matrix
#'
#' @return a data.frame with 1 columns containing the fold change


fc <- function(pair.mat){
  pair.mat[,2]/pair.mat[,1]
}


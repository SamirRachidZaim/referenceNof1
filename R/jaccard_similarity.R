#' random forest feature selection based on binomial exact test
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage jaccard(a,b)
#'
#' @param a a vector of DEG calls
#' @param b a vector of DEG calls
#'
#' @return a scalar indicating the jaccard similarity


jaccard <- function(a,b){
  num = length(intersect(a,b))
  denom= length(union(a,b))
  return(num/denom)
}
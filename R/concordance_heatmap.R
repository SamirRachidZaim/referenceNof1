#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage concordance_heatmap(jaccard_matrix)
#'
#' @param jaccard_matrix the concordance matrix used to create the heatmap
#'
#' @export
#'
#'
concordance_heatmap <- function(jaccard_matrix){
  mat <- jaccard_matrix
  mat2 <- data.frame(data.table::melt(mat))
  mat2$JaccardIndex = arules::discretize(mat2$value, method = 'fixed', 
                                 breaks = c(  0,0.5, 0.75, .9, 1))
  
  JaccardIndex <- mat2$JaccardIndex
  
  p <- ggplot2::ggplot(mat2, ggplot2::aes(mat2[,'Var1'], mat2[,'Var2'])) + 
    ggplot2::geom_raster(ggplot2::aes( fill=JaccardIndex)) 

  
  p <- p + ggplot2::theme_dark() + ggplot2::theme(
    axis.text.y = ggplot2::element_text(angle = 0, hjust = 1,size = ggplot2::rel(2)),
    title = ggplot2::element_text(angle = 0, hjust = 1,size = ggplot2::rel(1.5))) +
    ggplot2::ggtitle("Agreement After Expression Filter") + ggplot2::theme_light()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,size = ggplot2::rel(2)),
                   axis.text.y = ggplot2::element_text(angle = 0, hjust = 1,size = ggplot2::rel(2)),
                   title = ggplot2::element_text(angle = 0, hjust = .5,size = ggplot2::rel(1.5))) +
    ggplot2::xlab('')+
    ggplot2::ylab('') +    ggplot2::scale_fill_manual(values = c("grey", "yellow", "blue", "green"))
  


  print(p)
}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=T, message=FALSE, warning=FALSE-------------------------------------
require(locfit)
mcf7 = referenceNof1::mcf7
res <- referenceNof1::applyFilter(mcf7, expressionCutoff = 30, fold_change_window = c(0,.3), 1)
jaccard_matrix <- res$JaccardMatrix


## ----echo=T, message=FALSE, warning=FALSE,message=FALSE-----------------------
referenceNof1::concordance_heatmap(jaccard_matrix = jaccard_matrix)


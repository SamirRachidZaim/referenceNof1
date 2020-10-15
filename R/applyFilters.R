#' Constructing robust reference standards for Nof1 studies for precision medicine
#'
#' \code{referenceNof1} is the R implementation of the reference biomarker algorithm by (Zaim 2020)
#'
#' @usage applyFilter(gs.mat = utils::data("mcf7"), expressionCutoff =30,
#'                    fold_change_window=c(1, 1.5), FDR.cutoff=0.1)
#'
#' @param gs.mat a matrix of paired samples used to develop the reference standard.
#'               All paired samples must be adjacent
#' @param expressionCutoff a scalar positive integer > 1 indicating the minimum count value of RNA-seq
#' @param fold_change_window a tuple indicating the fold change range to search
#' @param FDR.cutoff fdr.threshold for determining which set of genes are DEGs in the reference standard
#'
#' @export


applyFilter <- function(gs.mat=utils::data("mcf7"), expressionCutoff=30, fold_change_window=c(1,1.5),FDR.cutoff=0.1){

  
  ### re-arrange
  ncols= ncol(gs.mat)
  gs.mat <- gs.mat[, c(seq(1, ncols,2),seq(2,ncols,2))]

  idx.cutoff <- which(rowMeans(gs.mat) > expressionCutoff)

  #### Filter FCs < 1.5 (i.e., effect size )
  #### before calculating ref standard degs

  gs.mat <- gs.mat[idx.cutoff, ]
  
  if(nrow(gs.mat)<100){
    stop("Error: less than 100 observations in selected cutoff")
  }

  #### Filter FCs < 1.5 (i.e., effect size )
  #### before calculating ref standard degs

  calculateFCperPair <- function(pair.mat){

    fc.mat <- data.table::data.table(FC = fc(pair.mat),
                         LogFC=log_fc(pair.mat))
    return(fc.mat)
  }

  ncols = ncol(gs.mat)

  fc.mat <- do.call(cbind,lapply(1:(ncols/2), function(x) calculateFCperPair(gs.mat[c(x,x+(ncols/2))])))
  fold_change.mat <- fc.mat[,seq(from = 1, to=ncols, by=2), with=F]
  log_fold_change.mat <- fc.mat[,seq(from = 2, to=ncols, by=2), with=F]

  fc.mean <- rowMeans(fold_change.mat)
  logfc.mean <- rowMeans(log_fold_change.mat)

  idx_between <- which(dplyr::between(abs(fc.mean), left = fold_change_window[1], right = fold_change_window[2]))

  #### Filter FCs < 1.5 (i.e., effect size )
  #### before calculating ref standard degs

  combine_both_idx = intersect(idx_between, idx.cutoff)

  gs.mat <- gs.mat[combine_both_idx, ]

  ############ ############ ############ ############ ############
  ############ ############ ############ ############ ############
  ############ ############ ############ ############ ############

  ## edgeR as Gold standard
  DEG_gs.edgeR <- .DEG_edgeR(countTable = gs.mat,
                            conditions =  c(rep('normal',(ncols/2)), rep('tumor',(ncols/2))),
                            replicateType = 'genetically identical model organisms',
                            disper = .1)

  y_true.edgeR <- p.adjust(DEG_gs.edgeR[[1]]$PValue, method='BY') < FDR.cutoff



  ## DESeq as Gold Standard
  gs.cohort.conditions =   c(rep('normal',(ncols/2)), rep('tumor',(ncols/2)))

  DESeq.preds  <- .DEG_DESeq.cb(gs.mat, conditions =gs.cohort.conditions )
  DESeq.preds <- DESeq.preds$padj
  DESeq.preds[which(is.na(DESeq.preds))] <- 1
  y_true.deseq <- DESeq.preds < FDR.cutoff


  ## DESeq2 as Gold Standard
  gs.cohort.conditions =   c(rep('normal',(ncols/2)), rep('tumor',(ncols/2)))

  DESeq2.preds  <- .DEG_DESeq2.cb(gs.mat, conditions =gs.cohort.conditions )
  DESeq2.preds <- DESeq2.preds$padj
  DESeq2.preds[which(is.na(DESeq2.preds))] <- 1

  y_true.deseq2 <- DESeq2.preds < FDR.cutoff

  ## DEGseq as Gold Standard
  gs.cohort.conditions =   c(rep('normal',(ncols/2)), rep('tumor',(ncols/2)))

  ## DEGseq
  DEGseq.preds <- .DEG_DEGseq.cb(gs.mat[,1:(ncols/2)],gs.mat[,((ncols/2)+1):ncols], FDR_threshold= .1)
  DEGseq.preds <- DEGseq.preds[order(DEGseq.preds$GeneNames),]
  DEGseq.preds <- DEGseq.preds[,'q.value.Benjamini.et.al..1995.']
  DEGseq.preds[which(is.na(DEGseq.preds))] <- 1
  y_true.degseq <- DEGseq.preds < FDR.cutoff #quantile(DEGseq.preds, .15) = 3.56227e-12

  #### NOISeq  goldstandard
  # countTable = ts.mat
  ## DESeq as Gold Standard
  gs.cohort.conditions =   c(rep('normal',(ncols/2)), rep('tumor',(ncols/2)))

  ## NOISeq
  NOISeq.cb.preds <- .DEG_noiseq.cb(gs.mat, conditions =gs.cohort.conditions)
  NOISeq.preds <- NOISeq.cb.preds$prob
  NOISeq.preds <-  p.adjust(1-NOISeq.preds, 'BY') #(transform to p-values)
  NOISeq.preds[which(is.na(NOISeq.preds))] <- 1
  y_true.noiseq <- NOISeq.preds < FDR.cutoff

  idx.nsq.idx <- which(y_true.edgeR ==T)
  idx.edg.idx <- which(y_true.deseq  ==T)
  idx.desq.idx <- which(y_true.deseq2 ==T)
  idx.desq2.idx <- which(y_true.degseq ==T)
  idx.degsq <- which(y_true.noiseq  ==T)

  y_true.intersection <- y_true.union <- logical(length(y_true.edgeR))

  y_true.intersection[Reduce(intersect, list(idx.nsq.idx, idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))] <- T
  y_true.union[Reduce(union, list(idx.nsq.idx, idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))] <- T


  lnint <- function(a,b){
    length(intersect(a,b))
  }

  lnint_pairwise <- function(a, b.list){
    sapply(b.list, function(x) lnint(a, x))
  }



  jaccard_pairwise <- function(a, b.list){
    sapply(b.list, function(x) jaccard_similarity(a, x))
  }

  nsq = lnint_pairwise(idx.nsq.idx, list(idx.nsq.idx, idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))
  edg = lnint_pairwise(idx.edg.idx, list( idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))
  dsq = lnint_pairwise(idx.desq.idx, list( idx.desq.idx,idx.desq2.idx,idx.degsq))
  dq2 = lnint_pairwise(idx.desq2.idx, list( idx.desq2.idx,idx.degsq))

  int.matrix = data.frame(NOISeq = numeric(5),
                          edgeR = numeric(5),
                          DESeq = numeric(5),
                          DESeq2= numeric(5),
                          DEGseq= numeric(5),
                          row.names = c('NOISeq','edgeR','DESeq','DESeq2','DEGseq'))

  int.matrix$NOISeq <- nsq
  int.matrix$edgeR[2:5] <- edg
  int.matrix$DESeq[3:5] <- dsq
  int.matrix$DESeq2[4:5] <- dq2
  int.matrix$DEGseq[5] <- length(idx.degsq)

  ## jaccard
  jac.nsq = jaccard_pairwise(idx.nsq.idx, list(idx.nsq.idx, idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))
  jac.edg = jaccard_pairwise(idx.edg.idx, list( idx.edg.idx,idx.desq.idx,idx.desq2.idx,idx.degsq))
  jac.dsq = jaccard_pairwise(idx.desq.idx, list( idx.desq.idx,idx.desq2.idx,idx.degsq))
  jac.dq2 = jaccard_pairwise(idx.desq2.idx, list( idx.desq2.idx,idx.degsq))



  jac.matrix = data.frame(NOISeq = numeric(5),
                          edgeR = numeric(5),
                          DESeq = numeric(5),
                          DESeq2= numeric(5),
                          DEGseq= numeric(5),
                          row.names = c('NOISeq','edgeR','DESeq','DESeq2','DEGseq'))

  jac.matrix$NOISeq <- jac.nsq
  jac.matrix$edgeR[2:5] <- jac.edg
  jac.matrix$DESeq[3:5] <- jac.dsq
  jac.matrix$DESeq2[4:5]<- jac.dq2
  jac.matrix$DEGseq[5] <- length(idx.degsq)

  jac.matrix <- jac.matrix + t(jac.matrix)
  diag(jac.matrix) <- rep(1,ncol(jac.matrix))

  int.matrix <- int.matrix + t(int.matrix)
  diag(int.matrix) <- rep(1,ncol(int.matrix))

  jac.matrix <- as.matrix(jac.matrix)

  return(list(JaccardMatrix = jac.matrix,
              Medians = robustbase::colMedians(jac.matrix))
  )

}

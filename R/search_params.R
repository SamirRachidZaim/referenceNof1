

# sample cutoffs
#  expressionCutoffs <- c(0, 30, 50)
#  FCRegions <- list(c(0,0.1), c(0.1, 0.2), c(0.2, 0.3),c(0.3, 0.5), c(0.5,Inf))


search_parameters <- function(expressionCutoffs, FCRegions,gs.mat, ts.mat, FDR.cutoff=0.1){

  params = expand.grid(expressionCutoffs, FCRegions)

  register()
  search_space <- mclapply(1:nrow(params),
                     function(x) applyFilter(expressionCutoff =  params[x,1],
                                                            fold_change_window=unlist(params[x,2]),
                                                            gs.mat=gs.mat,
                                                            ts.mat=ts.mat,
                                                            FDR.cutoff =FDR.cutoff))



}

##########################################################################################################################################################################################################################################################################################
##README: LMM-CAR function
##########################################################################################################################################################################################################################################################################################

#' @keywords internal
# Use CAR-LMM model for deconvolution
Decov <- function(ncpu,bulk,ExpCell,ExpBin,start,covariates,maxiter = maxiter,phi,CAR,cell_trajectory){

  ### construct the cell-state dependent covariance matrix (used in CAR)
  if(CAR){
    message('\n',"Computing the CAR matrix")
    ### register the environment
    BPPARAM = BiocParallel::MulticoreParam(workers = ncpu,progressbar = T)

    ### compute the distance of cells over the cell-state trajectory
    ED = rdist::rdist(as.matrix(cell_trajectory))
    
    ### gaussian-kernel transformation  
    kernel_mat =  as.matrix(exp(-ED^2 / (2)))
    D = diag(rowSums((kernel_mat)))
    CAR_mat = function(p){solve(D-p*kernel_mat)}
    S = BiocParallel::bplapply(FUN = CAR_mat,phi,BPPARAM = BPPARAM)
    
    ### store the calculared covariance matrix
    S[[length(S)+1]] = diag(ncol(ExpCell))
  }else{
    ### if not CAR, the covariance matrix will be the digonal matrix
    S = diag(ncol(ExpCell))
  }

  ### Run REML---(see /src/REML.cpp for the reml function) ------------------------------------------------------------------------------------------------
  ### Run fist bulk sample to determine the initial iteration value
  start = reml(start = start,X = as.matrix(covariates),y = as.matrix(bulk[,1]),Z = list(ExpCell) ,maxiter = maxiter,S=diag(nrow(cell_trajectory)))[[2]]
  message('\n',paste0(paste0('Run MeDuSA with ',ncpu)),' cores.')
  
  ### load the data to the environment
  ncpu = min(ncpu,parallel::detectCores())
  covariates = as.matrix(covariates)
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("bulk","reml","ExpCell","ExpBin","maxiter",'covariates',"start","S"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)

  ### show process
  pb = utils::txtProgressBar(min = 1, max = ncol(bulk), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  sampleID = NULL

  ### run LMM-CAR for each bulk sample
  abundance = foreach::foreach(sampleID = colnames(bulk), .options.snow = opts) %dopar2% {
    b = bulk[,sampleID]
    fixcmp = covariates
    rancmp = list(ExpCell)

    if(is.list(S)){
      ### reml-CAR iteration
      logL = -Inf
      for(i in seq_len(length(S))){
        s = S[[i]]
        parameter = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = rancmp ,maxiter = maxiter,S = s)
        start = parameter[[2]]
        ### chose the CAR-covariance matrix leading to the heighest log likely-hood  
        if(parameter[[3]]>logL){vi = parameter[[1]];logL = parameter[[3]]}
      }
    }else{
      ### only given the digonal CAR-covariance matrix (when CAR=FALSE)
      vi = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = rancmp ,maxiter = maxiter,S = S)[[1]]
    }

    ### compute coefficients of fixed terms using the generalized least squares
    beta = apply(ExpBin,2,function(x){
      x = cbind(x,fixcmp)
      (solve(t(x) %*% vi %*% x) %*% (t(x) %*% vi %*% b))[1]
    })
    return(beta)
  }

  ### log-out the environment
  parallel::stopCluster(cl)

  ### merge the estimated cell-state abundance
  abundance = do.call(cbind,abundance)


  return(abundance)
}

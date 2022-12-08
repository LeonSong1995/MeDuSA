#' @keywords internal
# Use CAR-LMM model for deconvolution
Decov = function(ncpu,bulk,g,ref,CBP,start,cov,maxiter=maxiter,phi,CAR,space){
  #1)---Build the covariance matrix------------
  if(CAR){
    BPPARAM = BiocParallel::MulticoreParam(workers=ncpu,progressbar = T)
    message('\n',"Computing the CAR matrix")
    ED = rdist::rdist(as.matrix(space))
    kernel_mat =  as.matrix(exp(-ED^2 / (2)))
    D = diag(rowSums((kernel_mat)))
    CAR_mat = function(p){solve(D-p*kernel_mat)}
    S = BiocParallel::bplapply(FUN = CAR_mat,phi,BPPARAM=BPPARAM)
    S[[length(S)+1]] = diag(ncol(ref))
  }else{
    S = diag(ncol(ref))
  }

  #2)---Run REML for the fist bulk sample to determine the initial iteration value
  start = reml(start = start,X = as.matrix(cov),y = as.matrix(bulk[g,1]),Z = list(as.matrix(ref[g,])) ,maxiter = maxiter,S=diag(nrow(space)))[[2]]

  #3)---Run deconvolution-----------------------
  message('\n',paste0(paste0('Run MeDuSA with ',ncpu)),' cores.')

  #3.1)--Load the data to the environment
  ncpu =ncpu
  cov = as.matrix(cov)
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("bulk","g","reml","ref","CBP","maxiter",'cov',"start","S"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(min = 1, max = ncol(bulk), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  sampleID = NULL

  #3.2)--Run CAR-LMM for each bulk sample
  abundance = foreach::foreach(sampleID = colnames(bulk), .options.snow = opts) %dopar2% {
    b = bulk[g,sampleID]
    fixcmp = cov
    rancmp = list(ref[g,])
    #3.2.1)--REML iteration (see REML.cpp for the reml function)
    if(is.list(S)){
      logL = -Inf
      for(i in seq_len(length(S))){
        s = S[[i]]
        parameter = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = rancmp ,maxiter = maxiter,S=s)
        start = parameter[[2]]
        if(parameter[[3]]>logL){vi = parameter[[1]];logL = parameter[[3]]}
      }
    }else{
      vi = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = rancmp ,maxiter = maxiter,S=S)[[1]]
    }
    #3.2.2)--Generalized least squares
    r = apply(CBP[g,],2,function(x){
      x=cbind(x,cov)
      (solve(t(x) %*% vi %*% x) %*% (t(x) %*% vi %*% b))[1]
    })
    return(r)
  }
  parallel::stopCluster(cl)
  abundance = do.call(cbind,abundance)
  return(abundance)
}

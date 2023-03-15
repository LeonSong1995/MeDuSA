#############################################################################################################
#' @title MeDuSA_VarExplain: the function to estimate explained variances of bulk RNA-seq data by the reference scRNA-seq data.

#' @param MeDuSA_obj The MeDuSA object. 
#' @param maxiter The iteration number of REML. Default by 1e+4.
#' @param ncpu The number of cpu cores to be used.
#'
#' @return \code{MeDuSA} returns: \itemize{
#' \item\code{MeDuSA_Object@VarianceExplain$VarExplain}: A numeric vector of the explained variances of bulk RNA-seq data by the reference scRNA-seq data.
#' \item\code{MeDuSA_Object@VarianceExplain$Chi_VarExplain}: A numeric vector of the corresponding squared-chi-value of the explained variances.
#' \item\code{MeDuSA_Object@VarianceExplain$Pval_VarExplain}: A numeric vector of the corresponding p-value of the explained variances.
#' }
#' @export

MeDuSA_VarExplain <- function(MeDuSA_obj,maxiter=1000,ncpu=4){

  ### check the input object
  if (!class(MeDuSA_obj) %in% "MeDuSA_Object"){
    stop("Please input the MeDuSA Object.")
  }

  ### scale the input matrix
  Y = scale(MeDuSA_obj@ExpProfile_bulk,center = T,scale = T)
  rancmp = scale(MeDuSA_obj@ExpProfile_cell,center = T,scale = T)
  fixcmp = MeDuSA_obj@covariates

  ### To accurately estimate the p-value of the explained variance, please do not use the CAR mode.
  S = as.matrix(diag(ncol(rancmp)))

  ### Run fist bulk sample to determine the initial iteration value
  start = reml(start = c(1e-3,1e-3),X = as.matrix(fixcmp),
               y = as.matrix(Y[,1]),Z = list(rancmp),maxiter = maxiter,S=S)[[2]]

  ### load the data to the environment
  ncpu = min(ncpu,parallel::detectCores())
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("Y","rancmp","fixcmp","S","reml","start"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)

  ### show process
  pb = utils::txtProgressBar(min = 1, max = ncol(Y), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  sampleID = NULL

  ### estimate explained variance and the corresponding p-value (chi-value)
  message('\n',paste0(paste0('Estimate the explained variances with ',ncpu)),' cores.')

  VarExp_Chi_Pvalue = foreach::foreach(sampleID = colnames(Y), .options.snow = opts) %dopar2% {
    b = Y[,sampleID]
    parameter = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = list(rancmp) ,maxiter = 1000,S = S)
    vi = parameter[[1]]
    varcmp = parameter[[2]]

    A = rancmp %*% t(rancmp) / (ncol(rancmp))
    Q = vi - vi %*% fixcmp %*% solve(t(fixcmp) %*% vi %*% fixcmp) %*% t(fixcmp) %*% vi
    Hi = solve(sum(diag(Q %*% A %*% Q %*% A)))
    chi = (varcmp[1]^2)/(2*Hi)
    Pvalue = pchisq(chi,df=1,lower.tail = F)
    VarExp = varcmp[1]/sum(varcmp)

    return(c(VarExp,chi,Pvalue))
  }
  parallel::stopCluster(cl)
  VarExp_Chi_Pvalue = as.data.frame(do.call(rbind,VarExp_Chi_Pvalue))

  ### update the MeDuSA_obj
  MeDuSA_obj@VarianceExplain$VarExplain = VarExp_Chi_Pvalue$V1
  MeDuSA_obj@VarianceExplain$Chi_VarExplain = VarExp_Chi_Pvalue$V2
  MeDuSA_obj@VarianceExplain$Pval_VarExplain = VarExp_Chi_Pvalue$V3
  names(MeDuSA_obj@VarianceExplain$VarExplain) =
  names(MeDuSA_obj@VarianceExplain$Chi_VarExplain) =
  names(MeDuSA_obj@VarianceExplain$Pval_VarExplain) = colnames(Y)

  return(MeDuSA_obj)
}

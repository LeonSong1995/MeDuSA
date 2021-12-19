# MeDuSA: mixed model-based deconvolution of cell-state abundance.
# Author: Liyang Song <songliyang@westlake.edu.cn>
#############################################################################################################
#' @title MeDuSA: mixed model-based deconvolution of cell-state abundance.
#' @description \code{MeDuSA} is a fine-resolution cellular deconvolution method, with the aim to use reference scRNA-seq data to predict the cell-state abundance of bulk RNA-seq data.
#'
#' @param bulk A matrix of bulk RNA-seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
#' @param sce A \code{seurat} object of the  reference scRNA-seq data (see \code{\link{seurat}}). Meta-data of the seurat object (sce@meta.data) need include two columns named as cell_type and cell_trajectory.
#' @param select.ct A character variable of the target cell type.
#' @param ncpu The number of CPU(s) to be used.
#' @param smooth A Boolean variable to determine whether to smooth the predicted cell-state abundance or not. The default value is TRUE.
#' @param smoothMethod A character variable to determine the smoothing method. The default method is \code{\link{loess}}.
#' @param gene A character vector of signature genes. The default value is NULL. With default, \code{MeDuSA} selects signature genes using the generalized additive model (see \code{\link{mgcv}}).
#' @param nbins A numeric variable to determine the number of bins along the cell trajectory, which is used to ensure the selected genes were uniformly scattered along the given trajectory. The default value is 10.
#' @param resolution A numeric variable to determine the number of cell bins along the cell trajectory (i.e., the resolution of the deconvolution analysis). The default value is 50.
#' @param knots A numeric variable to specify the number of knots. The default value is 10.
#' @param start A numeric vector for the initial value of the REML. Default by c(1e-5,1e-2).
#' @param maxgene A numeric variable to determine the number of signature genes. The default value is 200.
#' @param family A character variable to specify the distribution of GAM. See \code{\link{family.mgcv}} for a full list of what is available.
#' @param gcov A matrix (vector) of covariates in the GAM. The default value is NULL.
#' @param Xc A matrix (vector) of fixed covariates in the linear mixed model (i.e., covariates for predicting cell-state abundance). The default value is NULL.
#' @param maxiter The iteration number of REML. Default by 1e+4.
#' @param adj A Boolean variable to determine whether to include covariates when predicting the cell-state abundance.
#'
#' @return \code{MeDuSA} returns: \itemize{
#' \item\code{abundance}: A numeric matrix of cell-state abundance. Each row corresponds to a cell bin and each column corresponds to a sample.
#' \item\code{gene}: A character vector of the selected signature genes.
#' \item\code{PesudoTimeCellbin}: A numeric vector of the pseudo-time for each cell bin.
#' }
#'
#' @author Liyang Song <songliyang@westlake.edu.cn>
#' @examples
#' ##Library the package
#' library(MeDuSA)
#'
#' ##Load the test data:
#' data(ref)
#' data(cellType)
#' data(cellTrajectory)
#' data(bulk)
#'
#' ##Build the 'Seurat' obejct:
#' sce = CreateSeuratObject(ref)
#' sce$cell_type = cellType
#' sce$cell_trajectory = rep(0,ncol(sce))
#' sce$cell_trajectory[rownames(Trajectory)]=Trajectory
#'
#' ##Run MeDuSA (run with 6 courses):
#' csab = MeDuSA(bulk=bulk,sce=sce,select.ct='Epithelium',ncpu=6)

MeDuSA = function(bulk,sce,select.ct,ncpu=1,smooth=TRUE,smoothMethod='loess',gene=NULL,nbins=10,resolution=50,knots=10,start=c(1e-5,1e-2),maxgene=200,family='gaussian',gcov=NULL,Xc=NULL,maxiter=1e+4,adj=FALSE){

	#Checking the format of the input parameters
	message("Thanks for using MeDuSA to perform cell-state abundance deconvolution analyses.")
	if(!("Seurat" %in% class(sce))){
	  stop('Please input Seurat Object.')
	}
	if(!'cell_type' %in% colnames(sce@meta.data)){
	  stop('Do you forget to input the cell_type? Please check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!'cell_trajectory' %in% colnames(sce@meta.data)){
	  stop('Do you forget to input the cell_trajectory? Please check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!is.null(select.ct)){
	  if(is.na(table(sce$cell_type %in% select.ct)['TRUE'])){
	    stop('Do you forget to specify the cell type? Please check the select.ct.')
	  }
	}else{
	  stop('Do you forget to specify the cell type? Please check the select.ct.')
	}

	#Prepare the reference
	index = which(sce$cell_type==select.ct)
	space = as.matrix(sce$cell_trajectory[index])
	space = as.matrix(space[order(space,decreasing = F),])
	ref = as.matrix(sce@assays$RNA@counts[,rownames(space)])
	commGene = intersect(rownames(bulk),rownames(sce))
	bulk = bulk[commGene,]
	ref = ref[commGene,]


	#Select the genes
	nbins = min(nbins,ncol(ref))
	if(is.null(gene)){
	  g_chi = geneSelect(exprsData = ref,space=space,bulk=bulk, maxgene=maxgene,nbins=nbins,cov=gcov,family=family,k=knots,ncpu=ncpu)
	  g = g_chi$g
	  chi = g_chi$chi

	}else{g=Reduce(intersect,list(gene,rownames(bulk),rownames(sce)))}


	#Prepare the incidence matrix of fixed covariates.
	Xc_input = Xc[intersect(commGene,rownames(Xc)),]
	Xc = as.matrix(rep(1,length(commGene)))
	if(!is.null(Xc_input)){Xc = cbind(Xc,Xc_input)}

	#Include expression profiles of other cells as fixed covariates.
	if(length(index)<ncol(sce)){
		refALL = as.matrix(sce@assays$RNA@counts[commGene,-index])
		Xc = cbind(Xc,rowMeans(refALL[commGene,]))
	}

	#adjust the bulk RNA-seq data by covariates
	if(adj==TRUE){
		Xc = as.matrix(Xc[,-1])
		rownames(Xc)=commGene
		bulk_adj = sapply(1:ncol(bulk),function(i){
			Xcov = rowMeans(as.matrix(Xc[g,]))
			coef = lm(bulk[g,i]~Xcov+0)$coefficients
			#normalize the covariates coefficients
      coef[coef<0] = 0
			coef = coef/sum(coef+1e-100)
			bulk[g,i] = bulk[g,i]-(as.matrix(Xc[g,]) %*% as.vector(coef))
		})
		rownames(bulk_adj) = g
		colnames(bulk_adj) = colnames(bulk)
		bulk = bulk_adj
	}

	resolution = min(resolution, ncol(ref))
	bin = cluster(space,nbins = resolution)
	names(bin) = rownames(space)
	CBP = t(aggregate(t(ref[g,names(bin)]),by=list(bin),FUN=mean)[,-1])
	bmed = aggregate(space,by=list(bin),FUN=median)[,-1]
	ref = list(as.matrix(ref[g,]))

	#Run deconvolution with the MLM
	abundance = Decov(ncpu = ncpu,bulk=bulk,g=g,ref=ref,CBP=CBP,start=start,maxiter=maxiter)

	##for the data not converged
	abundance[,is.na(colSums(abundance))]=0

	#Smooth
	if(smooth==TRUE){
		if(smoothMethod=='loess'){
			abundance =apply(abundance,2,function(ab){predict(stats::loess(ab~bmed))})
		}else{
			num = max(5,round(length(bmed)*0.2))
			neighbour=sapply(1:length(bmed),function(i){order(abs(bmed[i]-bmed),decreasing = F)[1:num]})
			abundance = apply(abundance,2,function(x){sapply(1:ncol(neighbour),function(i){
				mean(x[neighbour[,i]])
			})})
		}
	}

	abundance = as.data.frame(abundance)

	###change the scale
	abundance = abundance/max(abundance*5)
	colnames(abundance) = colnames(bulk)
	rownames(abundance) = bmed
	return(list('abundance'=abundance,'gene'=g,'PesudoTimeCellbin'=bmed))
}


#' @keywords internal
#Mixed model deconvolution
Decov = function(ncpu,bulk,g,ref,CBP,start,maxiter=1e+4){
  ncpu =ncpu
  message('\n',paste0(paste0('Run deconvolution with ',ncpu)),' cores.')
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("bulk","g","reml","ref","CBP","maxiter","start"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(min = 1, max = ncol(bulk), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  geneNumber = NULL

  abundance = foreach::foreach(geneNumber = colnames(bulk), .options.snow = opts) %dopar2% {
    b = bulk[g,geneNumber]
    fixcmp = rep(1,length(g))
    vi = reml(start = start,X = as.matrix(fixcmp),y = as.matrix(b),Z = ref,maxiter = maxiter)[[4]]
    r = apply(CBP,2,function(x){solve(t(x) %*% vi %*% x) %*% (t(x) %*% vi %*% b)})
  }
  parallel::stopCluster(cl)
  abundance = do.call(cbind,abundance)
  return(abundance)
}




# MeDuSA: mixed model-based deconvolution of cell-state abundance.
# Author: Liyang Song <songliyang@westlake.edu.cn>


#############################################################################################################
#' @title MeDuSA: mixed model-based deconvolution of cell-state abundance。
#' @description \code{MeDuSA} is a fine-resolution cellular deconvolution method, with the aim to use reference scRNA-seq data to predict cell abundance distributed along a cell-state trajectory in a bulk RNA-seq data. MeDuSA is well suitable for biological scenarios in which the underlying mechanisms are associated with continuous transitions of cell-states.
#'
#' @param bulk A matrix of bulk RNA-Seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
#' @param sce A \code{seurat} object of the single-cell RNA-Seq data (see \code{\link{seurat}}). Meta data of this seurat object need include 'cellType' and 'cellTrajectory'.
#' @param selectCellType A character variable of the target cell type.
#' @param ncpu The number of CPUs to be used.
#' @param smooth A Boolean variable to determine whether to smooth the predicted cell-state abundance along the cell trajectory or not. The default value is TRUE. With default, the predicted cell-state abundance is smoothed using the Loess or averaging with its neighbors.
#' @param gene A character vector of signature genes. The default value is NULL. With default, \code{MeDuSA} selects signature genes to distinguish cells within different cell-states using the generalized additive model (GAM, see \code{\link{mgcv}}).
#' @param nbins A numeric variable to determine the number of gene bins along the cell trajectory, which is used to ensure the selected genes were uniformly scattered along the given trajectory. The default value is 10.
#' @param resolution A numeric variable to determine the number of cell bins along the cell trajectory (i.e., the resolution of the deconvolution analysis). The default value is 50.
#' @param knots A numeric variable to specify the numbers of knots of GAM. The default value is 10.
#' @param maxgene A numeric variable to determine the maximum number of genes to be selected. The default value is 200.
#' @param adj A Boolean variable to determine whether to correct covariates in predicting cell-state abundance or not.
#' @param gcov A matrix (or vector) of covariates in the GAM (i.e., covariates for selecting signature genes). The default value is NULL.
#' @param Xc  A matrix (or vector) of fixed covariates in the linear mixed model (i.e., covariates for predicting cell-state abundance). The default value is NULL. With default, MeDuSA includes expression profiles of other cell types as fixed covariates.
#' @param maxiter The maximum iterations of AI-REML. The default value is 1e+4.
#' @param smoothMethod  A character variable to specify the smoothing approach, including “loess” and “average”.
#' @param family A character variable to specify the distribution of GAM. See \code{\link{family.mgcv}} for a full list of what is available.

#' @details MeDuSA is a fine-resolution cellular deconvolution method that aims to use reference scRNA-seq data to predict cell abundance distributed along a cell-state trajectory in a bulk RNA-seq data set,
#' where cell trajectory specifies each cell as a point in a pseudo-time vector, determining the pattern of a dynamic process experienced by the cells.
#' \code{MeDuSA} comprises three steps:\itemize{
#' \item Select signature genes: \code{MeDuSA} associates genes with the cell trajectory using the generalized additive model (GAM). The association strength is scored using the Wald chi-squared value. MeDuSA then selects genes with the highest chi-squared values as signature genes.
#' \item Predicate cell-state abundance: \code{MeDuSA} predicts cell-state abundance for each cell bin using the linear mixed model.
#' \item Smooth cell-state abundance: \code{MeDuSA} smooths the predicted cell-state abundance by performing loess regression or averaging cells-state abundance with its neighbors.
#' }


#' @return \code{MeDuSA} returns: \itemize{
#' \item\code{abundance}: A matrix of cell-state abundance. Each row corresponds to a cell bin and each column corresponds to a sample.
#' \item\code{gene}:A vector of the selected signature genes.
#' \item\code{PesudoTimeCellbin}:A vector of the pseudo-time for each cell bin.
#' }
#'
#' @author Liyang Song <songliyang@westlake.edu.cn>
#' @examples
#' ##Library the package
#' library(MeDuSA)
#'
#' ##Load the test data
#' data(ref)
#' data(cellType)
#' data(cellTrajectory)
#' data(bulk)
#'
#' ##Build the 'Seurat' obejct:
#' sce = CreateSeuratObject(ref)
#' sce$cellType = cellType
#' sce$cellTrajectory = rep(0,ncol(sce))
#' sce$cellTrajectory[rownames(Trajectory)]=Trajectory
#'
#' ##Run MeDuSA (with 6 CPU cores)
#' CellStateAbundance = MeDuSA(bulk=bulk,sce=sce,selectCellType='Epithelium',ncpu=6)

MeDuSA = function(bulk,sce,selectCellType,ncpu=1,smooth=TRUE,smoothMethod='loess',gene=NULL,nbins=10,resolution=50,knots=10,maxgene=200,family='gaussian',gcov=NULL,Xc=NULL,maxiter=1e+4,adj=FALSE){

	#Checking the format of the input parameters
	message("Thanks for using MeDuSA to perform cell-state abundance deconvolution analysis.")
	if(!("Seurat" %in% class(sce))){
	  stop('Please input Seurat Object.')
	}
	if(!'cellType' %in% colnames(sce@meta.data)){
	  stop('Do you forget to input the cellType? Please check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!'cellTrajectory' %in% colnames(sce@meta.data)){
	  stop('Do you forget to input the cellTrajectory? Please check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!is.null(selectCellType)){
	  if(is.na(table(sce$cellType %in% selectCellType)['TRUE'])){
	    stop('No cell types selected. Please check the selectCellType')
	  }
	}else{
	  stop('No cell types selected. Please check the selectCellType')
	}


	#Prepare the reference
	index = which(sce$cellType==selectCellType)
	space = as.matrix(sce$cellTrajectory[index])
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
	Xc_input = Xc[commGene,]
	Xc = as.matrix(rep(1,length(commGene)))
	if(!is.null(Xc_input)){Xc = cbind(Xc,Xc_input)}

	#Include expression profiles of other cells as fixed covariates.
	if(length(index)<ncol(sce)){
		refALL = as.matrix(sce@assays$RNA@counts[,-index])
		Xc = cbind(Xc,rowMeans(refALL[commGene,]))
	}

	#adjust the bulk data by covariates
	if(adj==TRUE){
		rownames(Xc)=commGene
		bulk_adj = sapply(1:ncol(bulk),function(i){residuals(lm(bulk[commGene,i]~Xc[commGene,]))})
		rownames(bulk_adj) = commGene
		colnames(bulk_adj) = colnames(bulk)
		bulk = bulk_adj
	}

	resolution = min(resolution, ncol(ref))
	bin = cluster(space,nbins = resolution)
	names(bin) = rownames(space)
	CBP = t(aggregate(t(ref[g,names(bin)]),by=list(bin),FUN=mean)[,-1])
	bmed = aggregate(space,by=list(bin),FUN=median)[,-1]
	ref = list(as.matrix(ref[g,]))

	#Run deconvolution with the mixed model
	abundance = Decov(ncpu = ncpu,bulk=bulk,g=g,ref=ref,CBP=CBP,maxiter=maxiter)

	#Smooth
	if(smooth==TRUE){
		if(smoothMethod=='loess'){
			abundance =apply(abundance,2,function(ab){predict(stats::loess(ab~bmed))})
		}else{
			num = 5
			neighbour=sapply(1:length(bmed),function(i){order(abs(bmed[i]-bmed),decreasing = F)[1:num]})
			abundance = apply(abundance,2,function(x){sapply(1:ncol(neighbour),function(i){ mean(x[neighbour[,i]])})})
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
Decov = function(ncpu,bulk,g,ref,CBP,maxiter=1e+4){
  ncpu =ncpu
  message('\n',paste0(paste0('Run deconvolution with ',ncpu)),' cores.')
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("bulk","g","reml","ref","CBP","maxiter"),
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
    vi = reml(start = rep(var(b),length(ref)+1),X = as.matrix(fixcmp),y = as.matrix(b),Z = ref,maxiter = maxiter)[[4]]
    r = apply(CBP,2,function(x){solve(t(x) %*% vi %*% x) %*% (t(x) %*% vi %*% b)})
  }
  parallel::stopCluster(cl)
  abundance = do.call(cbind,abundance)
  return(abundance)
}


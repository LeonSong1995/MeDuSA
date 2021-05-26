# Mixed model-based Cell Trajectory Deconvolution (MCTD)
# Author: Liyang Song <songliyang@westlake.edu.cn>


#############################################################################################################
#' @title MCTD: Mixed model-based Cell Trajectory Deconvolution.
#' @description This is the main function of the "MCTD" method. \code{MCTD} is used to predict cell-abundance along a predefined cell trajectory in the given bulk RNA-seq data. \code{MCTD} is well suitable for biological scenarios in which the underlying mechanisms are associated with continuous transitions of cellular states.
#' @param bulk  A matrix of bulk RNA-Seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
#' @param sce  A \code{Seurat} object of the single-cell RNA-Seq data (see \code{\link{Seurat}}). Meta data of the 'sce' need include 'cellType' and 'cellTrjaectory'.
#' @param selectCellType A character of the names of the target cell-type.
#' @param ncpu The number of CPU cores to be used.
#' @param smooth A Boolean variable to determine whether to smooth the cell-abundance along the cell-trajectory or not. The default value is \emph{TRUE}. With default, the predicted cell abundance will be smoothed along the cell trajectory via the Loess (see \code{\link{loess}}).
#' @param gene A character vector of the genes to use as signatures.The default value is \emph{NULL}. With default, \code{MCTD} will select genes via the generalized additive model (GAM,see \code{\link{mgcv}}).
#' @param nbins A numeric variable to determine the number of cell-bins used in the gene selection. The parameter \code{nbins} was used to ensure the selected genes were evenly distributed along the cell trajectory. The default value is \emph{10}.
#' @param resolution A numeric variable to determine the resolution of deconvolution. After setting \code{resolution}, cells will be grouped into \emph{"resolution"} cell-bins. The default value is \emph{50}.
#' @param knots A numeric variable to specify the numbers of knots used in the GAM. The default value is \emph{10}.
#' @param maxgene A numeric variable to specify the maximum numbers of genes to be selected. The default value is \emph{200}.
#' @param mode A character variable to specify the distribution and link to use in the gene selection. See \code{\link{family.mgcv}} for a full list of what is available. Note that we suggest to use the \emph{gaussian} or the \emph{nb} (negative binomial distribution)' here.The default mode is the \emph{gaussian}.
#' @param gcov A matrix (or vector) of fixed covariates in gene selection,such as donors, genders and extra.The default value is \emph{NULL}.
#' @param Xc  A matrix (or vector) of fixed covariates in deconvolution.The default value is \emph{NULL}. With default, \code{MCTD} will include expression profiles of other cell-types as the deconvolution covariates.
#' @param maxiter The maximum iterations of REML.The default value is \emph{1e+4}.
#'
#' @return \code{MCTD} returns: \itemize{
#' \item\code{abundance}: A matrix of cell-abundance. Each row corresponds to a median trajectory point of the cell-bins and each column corresponds to a certain sample.
#' \item\code{gene}: A vector of the signature genes used in deconvolution.
#' }
#'
#' @details  Unlike most deconvolution methods decomposing cell-abundance in bulk tissues for cell-types,
#' \code{MCTD} is a fine-resolution deconvolution method that predicts cell abundance along a predefined cell trajectory.
#' Cell-trajectory specifies each cell as a point in a geometric vector, determining the pattern of a dynamic process experienced by cells.
#' With inputting cell trajectory and scRNA-seq as references, \code{MCTD} infers the cell-abundance along the cell trajectory within the given bulk RNA-seq data.
#' \code{MCTD} contains three steps:\itemize{
#' \item Gene selection: \code{MCTD} selects signature genes to distinguish cells along the cell-trajectory via the GAM.
#' \item Deconvolution: \code{MCTD} infers the cell-abundance for each cell-bin via the mixed model.
#' \item Smoothing: \code{MCTD} smooths the cell-abundance along the cell-trajectory via the loess.This step is performed if smooth = TRUE.}
#'
#' @author Liyang Song <songliyang@westlake.edu.cn>
#' @examples
#' ##Load the test data
#' data(ref)
#' data(cellType)
#' data(cellTrjaectory)
#' data(bulk)
#'
#' ##Build the 'Seurat' obejct:
#' sce = CreateSeuratObject(ref)
#' sce$cellType = cellType
#' sce$cellTrjaectory = rep(0,ncol(sce))
#' sce$cellTrjaectory[rownames(Trajectory)]=Trajectory
#'
#' ##Run MCTD (with 6 CPU cores)
#' CellAbundance = MCTD(bulk=bulk,sce=sce,selectCellType='Epithelium',ncpu=6)$abundance

MCTD = function(bulk,sce,selectCellType,ncpu=1,smooth=TRUE,gene=NULL,nbins=10,resolution=50,knots=10,maxgene=200,mode='gaussian',gcov=NULL,Xc=NULL,maxiter=1e+4){

	#Checking the format of the input parameters
	message("Thanks for using MTD to perform cell-trjaectory deconvolution analysis.")
	if(!("Seurat" %in% class(sce))){
	  stop('Please input Seurat Object.')
	}
	if(!'cellType' %in% colnames(sce@meta.data)){
	  stop('Please input cellType, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!'cellTrjaectory' %in% colnames(sce@meta.data)){
	  stop('Please input cellTrjaectory, check your MetaData: Seurat_Obj@meta.data.')
	}
	if(!is.null(selectCellType)){
	  if(is.na(table(sce$cellType %in% selectCellType)['TRUE'])){
	    stop('No cell-types selected. Please check the selectCellType')
	  }
	}else{
	  stop('No cell-types selected. Please check the selectCellType')
	}


	#Prepare the reference
	index = which(sce$cellType==selectCellType)
	space = as.matrix(sce$cellTrjaectory[index])
	space = as.matrix(space[order(space,decreasing = F),])
	Ref = as.matrix(sce@assays$RNA@counts[,rownames(space)])
	bulk = bulk

	#Select the genes
	nbins = min(nbins,ncol(Ref))
	if(is.null(gene)){
	  g = geneSelect(exprsData = Ref,space=space,bulk=bulk, maxgene=maxgene,nbins=nbins,cov=gcov,mode=mode,k=knots,ncpu=ncpu)$g
	}else{g=gene}

	#Prepare the incidence matrix of fixed covariates.
	Xc_input = Xc[g,]
	Xc = as.matrix(rep(1,length(g)))
	if(!is.null(Xc_input)){Xc = cbind(Xc,Xc_input)}

	#Include expression profiles of other cell-types as fixed covariates.
	if(length(index)<ncol(sce)){
	  temp = sce[g,-index]
	  Xc = cbind(Xc,t(aggregate(t(as.matrix(temp@assays$RNA@counts)),by=list(temp$cellType),FUN=mean)[,-1]))
	}
	rownames(Xc)=g

	#Prepare the cell-trajectory bin and cell-bin expression profiles
	resolution = min(resolution, ncol(Ref))
	bin = cluster(space,nbins = resolution)
	names(bin) = rownames(space)
	CBP = t(aggregate(t(ref[g,names(bin)]),by=list(bin),FUN=mean)[,-1])
	bmed = aggregate(space,by=list(bin),FUN=median)[,-1]

	#Run deconvolution with the mixed model
	abundance = Decov(ncpu = ncpu,bulk=bulk,g=g,cov=Xc,ref=ref,CBP=CBP,maxiter=maxiter)

	#Smooth
	if(smooth==TRUE){
	  abundance =apply(abundance,2,function(ab){predict(stats::loess(ab~bmed))})
	}
	abundance = as.data.frame(abundance)
	colnames(abundance) = colnames(bulk)
	rownames(abundance) = bmed
	return(list('abundance'=abundance,'gene'=g))
}


#' @keywords internal
#Mixed model deconvolution
Decov = function(ncpu,bulk,g,cov,ref,CBP,maxiter=1e+4){
  ncpu =ncpu
  message('\n',paste0(paste0('Run deconvolution with ',ncpu)),' cores.')
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("bulk","g","reml","cov","ref","CBP","maxiter"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(min = 1, max = ncol(bulk), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  geneNumber = NULL

  abundance = foreach::foreach(geneNumber = colnames(bulk), .options.snow = opts) %dopar2% {
    b = bulk[g,geneNumber]
    fixcmp = cov[g,]
    vi = reml(start = c(1,1),X = as.matrix(fixcmp),y = as.matrix(b),Z = list(as.matrix(ref[g,])),maxiter = maxiter)[[4]]
    r = apply(CBP,2,function(x){solve(t(x) %*% vi %*% x) %*% (t(x) %*% vi %*% b)})
  }
  parallel::stopCluster(cl)
  abundance = do.call(cbind,abundance)
  return(abundance)
}


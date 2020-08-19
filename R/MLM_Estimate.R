# Functions to perform bulk RNA-seq deconvolution analysis with the linear mixed model
#
# Author: Liyang Song <liyang.song@ifar.ac.cn>
# Advisor: Jian Yang, Xiwei Sun
# Copyright: Liyang Song
#############################################################################################################
#' Bulk RNA-seq deconvolution analysis with the linear mixed model
#'
#' @param bulk: matrix of bulk tissue expression;
#' @param sce: object of SingleCellExperiment, single-cell reference data;
#' @param gene: vector of signature genes;
#' @param data_type: character, type of the inputed single-cell reference data (count/tpm/rpkm);
#' @param select.ct: vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x};
#' @param RanSplit: character, the label for how to separate random components. Default as \code{NULL}. If \code{NULL}, all cells will be fitted as one random component;
#' @param ct.cell.size: vector of cell sizes with labeled cell type names. Default as NULL. If NULL, then estimate cell size from data;
#' @param BatchCorrect: bool, whether to remove the batch effect between bulk data and single-cell reference data or not. Default as FALSE;
#' @param Filter: bool, whether to remove outlier cells in single-cell reference data or not. Default as TRUE;
#' @param SF: numeric, scaling factor.Default as 1e+3;
#' @param DoParallel: bool, whether to perform analysis parallelly or not. Default as FALSE;
#' @param ncpu: numeric, the number of CPUs used. Default as NULL. If NULL, then maximum detected CPUs - 1 will be used ;
#' @param verbose: bool, whether to show running information or not during deconvolution analysis. Default as FALSE;
#' @param iter_max: numeric, maximum iteration number.Default as 1e+3.
#'
#' @return a list with elements:
#'   *ct.pro: matrix of cell type proportions estimated by mixed linear model (sample x cell type);
#'   *ct.pro.p: matrix of p value (χ²(df=1)) for the cell type proportions estimated by the mixed linear model (sample x cell type);
#'   *cellSize: vector of cell sizes with labeled cell type names.
#'
#' @export
#'
#' @examples
#' library(MLM)
#' ct.es = Estimate(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type = 'count')
#' cor(ct.es$ct.pro[,colnames(ct.real)],ct.real)
#'
Estimate = function(bulk,sce,gene,data_type, select.ct = NULL, RanSplit=NULL, ct.cell.size = NULL,BatchCorrect=F,Filter=T,SF=1e+3,DoParallel=FALSE, ncpu=NULL, verbose=FALSE,iter_max=1000){

	message("Thanks for using MLM to perform bulk deconvolution analysis.\n")

	# Checking the correct format of the reference single-cell data input
	if (!(data_type %in% c('count','tpm','rpkm'))){
		stop('Please input correct data_type: cout/tpm/rpkm.')
	}
	if(!("SingleCellExperiment" %in% class(sce))){
		stop('Please input correct single-cell data: "SingleCellExperiment"(https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).')
	}
	if(!is.null(RanSplit)){
		if (is.na(match(RanSplit,colnames(SingleCellExperiment::colData(sce))))){stop('Do not know how to split randomp components, please check your MetaData: colData(sce).')}
	}
	if(!'sampleID' %in% colnames(SingleCellExperiment::colData(sce))){
		stop('Please input sampleID, check your MetaData: colData(sce).')
	}
	if(!'cellType' %in% colnames(SingleCellExperiment::colData(sce))){
		stop('Please input cellType, check your MetaData: colData(sce).')
	}
	if(length(unique(SingleCellExperiment::colData(sce)$cellType))==1){
		stop('MLM can not work with only one cell type inputted.')
	}
	if(!is.null(select.ct)){
		if(is.na(table(sce$cellType %in% select.ct)['TRUE'])){
			stop('No cell types selected. Please check the select.ct!')
		}else{
			sce = sce[,sce$cellType %in% select.ct]
		}
	}

	MetaData = SingleCellExperiment::colData(sce)
	exprsData = assay(sce,data_type)


	# Checking the signature genes input
	commonGene = intersect(rownames(exprsData),rownames(bulk))
	gene = intersect(commonGene,gene)
	if(length(gene)<10){
		stop('Too few signature genes (signature genes < 10).')
	}

	bulk = as.matrix(bulk)
	MetaData$cellType = as.vector(MetaData$cellType)
	MetaData$sampleID = as.vector(MetaData$sampleID)

	# Preparing for the basic running information (fixed/random components, cell size...)
	message("Data preparing...")
	Info = basis(bulk = bulk, exprsData = exprsData, MetaData=MetaData, ct.cell.size = ct.cell.size,
				data_type=data_type,gene=gene,BatchCorrect = BatchCorrect,Filter=Filter,SF=SF)


	base = Info$base
	bulk = Info$bulk
	data_cellType = Info$data_cellType
	cellSize = Info$cellSize
	type_n = length(data_cellType)

	# Preparing for the random components
	rancmp = sapply(seq(1,type_n), function(i){
				Z_new = data_cellType
				Z_new[i] = NULL
				Rest =  matrix(unlist(Z_new),nrow = length(gene))
				cellName = unlist(lapply(Z_new,colnames))
				colnames(Rest) = cellName
				Rest})

	ct_name = colnames(base)
	bk_name = colnames(bulk)

	# Checking the cpus
	if (DoParallel){
	  # File for output the warning information
		path = paste(getwd(),'RunLog',sep = '/')
		dir.create(path)
		filename = paste(path,paste(Sys.time(), paste(Sys.getpid(),'logger',sep='_'),sep='_'),sep='/')

		if(is.null(ncpu) | (detectCores()-1)<ncpu ){
			ncpu = detectCores()-1
		}
		message(paste0("Please wait, running bulk deconvolution with ",ncpu,' CPUs...'))
		cl = makeCluster(ncpu)
		registerDoParallel(cl)

		# Run REML parallelly
		estimate = parSapply(cl,1:ncol(bulk),function(sid){
					r = RunReml(sid,base=base,bulk=bulk,rancmp=rancmp,MetaData=MetaData, RanSplit=RanSplit, iter_max,ct_name, bk_name,verbose,filename,DoParallel)
					b = pmax(r['b',],0)
					p = r['p',]
					return(c(b,p))})
		stopCluster(cl)
		cat("Deconvolution finished.",file=filename,append=T)

	}else{
		message("Please wait, running bulk deconvolution...")
	  # Run REML
		estimate = sapply(1:ncol(bulk),function(sid){
					r = RunReml(sid,base=base,bulk=bulk,rancmp=rancmp,MetaData=MetaData, RanSplit=RanSplit, iter_max,ct_name, bk_name,verbose,filename,DoParallel)
					b = pmax(r['b',],0)
					p = r['p',]
					return(c(b,p))})
	}

	b = t(estimate[seq(1:type_n),])
	p = t(estimate[-seq(1:type_n),])
	colnames(b) = colnames(p) = ct_name
	rownames(b) = rownames(p) = bk_name

	# Adjusting for the cell size when count matrix inputted
	if(data_type=='count'){
		b = sweep(b,2,cellSize[colnames(b)],'/')
		b = sweep(b,1,rowSums(b),'/')
	}else{
		warning("The estimated cell type proportions are not comparable among different cell types!")
	}
	return(list('ct.pro'=b,'ct.pro.p'=p,'cellSize'= cellSize))

}




#' Run residual maximum likelihood (REML) to fit the linear mixed model.
#'
#' @param sid: numeric, id of the inputed samples;
#' @param base: matrix of the fixed components (gene x cell type);
#' @param bulk: matrix of the bulk RNA-seq samples (gene x sample);
#' @param rancmp: matrix of the random components (gene x cell);
#' @param MetaData: data frame of the summary information of bulk RNA-seq samples. Include sampleID, cellTyp;
#' @param RanSplit: character, label for how to seperate random components. Default as \code{NULL}. If \code{NULL}, all cells will be fitted as one random component;
#' @param iter_max: numeric, maximum iteration number.Default as 1e+3;
#' @param ct_name: vector of cell type names;
#' @param bk_name: vector of bulk RNA-seq sample names;
#' @param verbose: bool, whether to show running information or not during deconvolution analysis;
#' @param filename: path, save the warning information during REML analysis;
#' @param DoParallel:  bool, whether to perform analysis parallelly.
#'
#' @return a list with elements:
#'   *b: matrix of cell-type proportions estimated by mixed linear model (sample x cell type);
#'   *p: matrix of p value (χ²(df=1)) for the cell type proportions estimated by mixed linear model (sample x cell type).
#'
#' @export
#'
#'
RunReml = function(sid,base, bulk, rancmp, MetaData, RanSplit, iter_max,ct_name, bk_name, verbose,filename=NULL,DoParallel){
  x = as.matrix(base)
  y = bulk[,sid]
  type_n = ncol(x)

  result = sapply(seq(1,type_n),function(id){
    logger = paste(paste('sample',bk_name[sid],sep = '_'),ct_name[id],sep = '_')
    fixcmp = x[,id]
    if (!is.null(RanSplit)){
      # Splitting random components based on the inputted RanSplit
      sepInfo = MetaData[intersect(colnames(rancmp[[id]]),rownames(MetaData)),RanSplit]
      Z = sapply(unique(sepInfo),function(sid){
        rancmp[[id]][,sepInfo %in% sid]
      })
    }else{
      Z = list(rancmp[[id]])
    }
    mlmfit = reml(X = fixcmp,y = y,Z = Z,maxiter = iter_max)
    #reml(X,y,Z,maxiter) function was written by c++ (see ~/src/REML.cpp)
    #mlmfit[1] : effect size for the fixed component;
    #mlmfit[2] : χ²(df=1) value for the fixed component;
    #mlmfit[3] : flag for reml convergence (0:not converge, 1: converge);
    #mlmfit[4] : flag for V matrix (0:not invertible, 1: invertible);
    #mlmfit[5] : flag for  X^t * V^-1 * X matrix (0:not invertible, 1: invertible);
    #mlmfit[6] : flag for iteration (0:not reach to the maximum iteration , 1:reach to the maximum iteration).
	if (DoParallel){
		if(mlmfit[3]){

		  b = mlmfit[1]
		  P_value=1-pchisq(mlmfit[2],1)

		  if(!mlmfit[6]){

			cat(paste(logger,"Log-likelihood not converged. The results are not reliable. You can specify the parameter max_iter to allow for more iterations.\n",sep = ':'),file=filename,append=T)
		  }else if (verbose){
			cat(paste(logger,'REML converged and estimation finished.\n',sep = ':'),file=filename,append=T)
		  }
		}else{
		  lmfit = as.vector(summary(lm(y~base+0))$coefficients[id,c('Estimate','Pr(>|t|)')])
		  b = lmfit[1]
		  P_value = lmfit[2]

		  if(!mlmfit[4]){
			cat(paste(logger,"The V matrix is not positive. Switch to multiple linear regression!\n",sep = ':'),file=filename,append=T)
		  }
		  if(!mlmfit[5]){
			cat(paste(logger,"X^t * V^-1 * X matrix is not invertible. Switch to multiple linear regression!\n",sep = ':'),file=filename,append=T)
		  }
		}
	}else{
		if(mlmfit[3]){

		  b = mlmfit[1]
		  P_value=1-pchisq(mlmfit[2],1)

		  if(!mlmfit[6]){

			warning(paste(logger,"Log-likelihood not converged. The results are not reliable. You can specify the parameter max_iter to allow for more iterations..",sep = ':'))
		  }else if (verbose){
			message(paste(logger,'REML converged and estimation finished.',sep = ':'))
		  }
		}else{
		  lmfit = as.vector(summary(lm(y~base+0))$coefficients[id,c('Estimate','Pr(>|t|)')])
		  b = lmfit[1]
		  P_value = lmfit[2]

		  if(!mlmfit[4]){
			warning(paste(logger,"The V matrix is not positive. Switch to multiple linear regression!",sep = ':'))
		  }
		  if(!mlmfit[5]){
			warning(paste(logger,"X^t * V^-1 * X matrix is not invertible. Switch to multiple linear regression!",sep = ':'))
		  }
		}
	}

    return(c('b'=b,'p'=P_value))
  })

  colnames(result) = colnames(x)
  return(result)
}

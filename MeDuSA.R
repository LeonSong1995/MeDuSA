# MeDuSA: mixed model-based deconvolution of cell-state abundance.
# Author: Liyang Song <songliyang@westlake.edu.cn>
#############################################################################################################
#' @title MeDuSA: mixed model-based deconvolution of cell-state abundance.
#' @description \code{MeDuSA} is a cellular deconvolution method that leverages scRNA-seq data as a reference to estimate cell-state abundance in bulk RNA-seq data.
#'
#' @param bulk A matrix of bulk RNA-seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
#' @param sce A \code{seurat} object of the  reference scRNA-seq data (see \code{\link{seurat}}). Meta-data of the seurat object (sce@meta.data) need include two columns named as cell_type and cell_trajectory.
#' @param select.ct A character variable of the target cell type.
#' @param ncpu The number of cpu cores to be used.
#' @param smooth A Boolean variable to determine whether to smooth the predicted cell-state abundance or not. The default value is TRUE.
#' @param smoothMethod A character variable to determine the smoothing method. The default method is \code{\link{loess}}.
#' @param nbins A numeric variable to determine the number of bins along the cell trajectory, which is used to ensure the selected genes were uniformly scattered along the given trajectory. The default value is 10.
#' @param resolution A numeric variable to determine the number of cell bins along the cell trajectory (i.e., the resolution of the deconvolution analysis). The default value is 50.
#' @param knots A numeric variable to specify the number of knots. The default value is 10.
#' @param start A numeric vector for the initial value of the REML. Default by c(1e-5,1e-2).
#' @param family A character variable to specify the distribution of gam See \code{\link{family.mgcv}} for a full list of what is available.
#' @param maxiter The iteration number of REML. Default by 1e+4.
#' @param adj A Boolean variable to determine whether to include covariates when predicting the cell-state abundance.
#' @param CAR A Boolean variable to determine whether to model abundance correlations among cells.
#' @param phi A numeric vector for searching the optimal cell correlations. The default value is c(0.2,0.4,0.6,0.9)
#' @param fixCov A matrix (vector) of fixed covariates in the model (i.e., covariates for estimating cell-state abundance). The default value is NULL.
#' @param markerGene A character vector of marker genes over the cell-state trajectory. The default value is NULL. With default, \code{MeDuSA} selects  genes using the gam-wald or wilcox test.
#' @param method A character variable to specify the method used in selecting marker genes (wilcox or gam).Default by wilcox.
#' @param geneNumber A numeric variable to determine the number of signature genes. The default value is 200.
#' @param span  A numeric variable to control the degree of smoothing \code{\link{loess}}.
#' @param neighbor A numeric variable to determine the number of neighboring cells used in smoothing (do not used when the smooth method is loess).
#' @param fractional A Boolean variable to determine whether to normalize the estimated cell-state abundance to the fractional abundance (0-1).Default by FALSE.
#'
#' @return \code{MeDuSA} returns: \itemize{
#' \item\code{MeDuSA_Object}: The MeDuSA object.
#' \item\code{MeDuSA_Object@Estimation$cell_state_abundance}: A matrix of cell-state abundance. Each row corresponds to a cell-state bin and each column corresponds to a bulk sample.
#' \item\code{MeDuSA_Object@Estimation$markerGene}: A character vector of the used marker genes.
#' \item\code{MeDuSA_Object@Estimation$TimeBin}: A numeric vector of the median pseudo-time for each cell-state bin.
#' }
#'
#' @author Liyang Song <songliyang@westlake.edu.cn>
#' @examples
#' see: https://github.com/LeonSong1995/MeDuSA

MeDuSA <- function(bulk,sce,select.ct=NULL,resolution=50,fixCov=NULL,adj=FALSE,
                   markerGene=NULL,nbins=10,knots=10,method="wilcox",family='gaussian',geneNumber=200,
                   CAR=FALSE,phi=c(0.2,0.4,0.6,0.9),
                   ncpu=1,start=c(1e-5,1e-2),maxiter=1e+4,
                   smooth=TRUE,smoothMethod='loess',span=0.35,neighbor=5,fractional=FALSE){

  #1)---Check inputs
  stopifnot(class(sce) == "Seurat")
  check_metadata(sce@meta.data, "cell_type")
  check_metadata(sce@meta.data, "cell_trajectory")

  if (is.null(select.ct)) {
    stop("Please specify the cell type using the 'select.ct' argument.")
  } else if (sum(sce$cell_type %in% select.ct) == 0) {
    stop("No cells were found with the specified cell type.")
  }
  message('MeDuSA: mixed model-based deconvolution of cell-state abundance')

  #2)---Select (Check) marker genes
  focal_cell_index = which(sce$cell_type == select.ct)
  if(is.null(markerGene)) {
    message(paste0("\n",paste0('Marker genes are not provided. MeDuSA will select marker genes over the cell trajectory.')))
    markerGene = MeDuSA_marker(sce[,focal_cell_index],bulk,
                               geneNumber = geneNumber,nbins = nbins,
                               family = family,k = knots,ncpu = ncpu,method = method)
  }else{
    markerGene = Reduce(intersect,list(markerGene,rownames(bulk),rownames(sce)))
  }

  #3)---Prepare fixed covariates of other cell types
  covariates = prepare_fixCov(sce = sce, markerGene = markerGene,fixCov = fixCov,
                              adj = adj,focal_cell_index = focal_cell_index)


  #4)---Prepare input matrix of the focal cell type
  focal_Ct = prepare_focal(sce,resolution,markerGene,focal_cell_index)
  TimeBin = focal_Ct$TimeBin
  ExpBin = focal_Ct$ExpBin
  ExpCell = focal_Ct$ExpCell
  cell_trajectory = focal_Ct$cell_trajectory
  remove(focal_Ct)

  #5)---Run LMM-CAR
  timeStart = Sys.time()
  abundance = Decov(ncpu = ncpu,
                    bulk = bulk[markerGene,],
                    ExpBin = as.matrix(ExpBin[markerGene,]),
                    ExpCell = as.matrix(ExpCell[markerGene,]),
                    cell_trajectory = as.matrix(cell_trajectory),
                    covariates = as.matrix(covariates[markerGene,]),
                    start = c(1e-3,1e-3),maxiter = maxiter,phi = phi,CAR = CAR)
  timeEnd = Sys.time()
  timeCost = difftime(timeEnd, timeStart, units='mins')
  message(paste("\n",paste(paste0('Elapsed time of LMM-CAR for ',paste0(ncol(bulk),' bulk-samples')),paste0(timeCost,' mins'),sep=':')))

  #6)---Set the abundance of not converging sample as 0
  abundance[,is.na(colSums(abundance))] = 0

  #7)---Smooth the estimated cell state abundance (optional)
  if(smooth){
    abundance = smooth_abundance(abundance = abundance,TimeBin = TimeBin,
                                 smoothMethod = smoothMethod,neighbor = neighbor,span = span)
  }

  #8)---Normalize the abundance to fractional abundance [0,1] (optional)
  if(fractional){
    abundance = sweep(abundance,2,abs(apply(abundance,2,min)),'+')
    abundance = sweep(abundance,2,colSums(abundance),'/')
  }

  #9)---Return (the MeDuSA object)
  colnames(abundance) = colnames(bulk)
  rownames(abundance) = colnames(ExpBin)

  MeDuSA_obj = Make_MeDuSA_Object(abundance = abundance,
                           TimeBin = TimeBin,
                           markerGene = markerGene,
                           ExpCell = ExpCell,
                           ExpBulk = bulk[markerGene,],
                           covariates = covariates)


  return(MeDuSA_obj)
}

##########################################################################################################################################################################################################################################################################################
##README: Base function of MeDuSA
##########################################################################################################################################################################################################################################################################################
### Includes:
### ----------------------------------------------------------------------------------------------------------------------
### 1) define and create the MeDuSA object. 
### 2) check_metadata: check the metadata of the given Seurat object (cell_type and cell_trajectory are requried).
### 3) prepare_fixCov: Prepare the matrix of fixed covariates.
### 4) prepare_focal: Prepare the matrices of focal cell type, which will return: 
###    1.median state-score(pseudo-time) cell state bin
###    2.mean expression profile of cell state bin
###    3.expression profile of individual cells
###    4.cell trajectory
### 5) smooth_abundance: Smooth the estimated cell-state abundance.
### ----------------------------------------------------------------------------------------------------------------------

#' @keywords internal
#' define the MeDuSA object
setClass("MeDuSA_Object",
         representation(
           Estimation = "list",
           ExpProfile_cell = "matrix",
           ExpProfile_bulk = 'matrix',
           covariates = "matrix",
           VarianceExplain = "list")
)

#' @keywords internal
#' create the MeDuSA object
Make_MeDuSA_Object <- function(abundance = NULL,TimeBin = NULL,markerGene = NULL,
                            ExpCell = NULL,ExpBulk = NULL,covariates = NULL,
                            VarExplain = NULL,Chi_VarExplain = NULL, Pval_VarExplain = NULL)
{

  MeDuSA_obj = new("MeDuSA_Object",
            Estimation = list('cell_state_abundance'= abundance,
                              'TimeBin'= TimeBin,
                              'markerGene'= markerGene),
            ExpProfile_cell = as.matrix(ExpCell),
            ExpProfile_bulk = as.matrix(ExpBulk),
            covariates = as.matrix(covariates),
            VarianceExplain = list('VarExplain' = VarExplain,
                                   'Chi_VarExplain' = Chi_VarExplain,
                                   'Pval_VarExplain' = Pval_VarExplain))
  return(MeDuSA_obj)
}

#' @keywords internal
#' check the metadata of the Seurat object
check_metadata <- function(metadata, column){
  if (!(column %in% colnames(metadata))) {
    stop(paste0("Missing '", column, "' column in @meta.data."))
  }
}

#' @keywords internal
#' Prepare the matrix of fixed covariates
prepare_fixCov <- function(sce,markerGene,fixCov,adj,focal_cell_index){
  if(adj==TRUE){
    ### Condition-1: Compute the fixed covariates using other cell types
    if(is.null(fixCov)){
      if(ncol(sce) > length(focal_cell_index)){
        fixCov = Seurat::AverageExpression(sce[markerGene,-focal_cell_index],group.by = 'cell_type',slot = 'counts',assays = 'RNA')[[1]]
        fixCov = cbind(1,fixCov)
      }else{
        fixCov = as.matrix(rep(1,length(markerGene)))
      }
    }else{
    ### Condition-2: Using fixed covariates given by users
      fixCov = as.matrix(cbind(rep(1,length(markerGene)),fixCov[markerGene,,drop=F]))
    }
  }else{
    ### Condition-3: No fixed covariates
    fixCov = as.matrix(rep(1,length(markerGene)))
  }

  fixCov = as.matrix(fixCov)
  rownames(fixCov) = markerGene
  return(fixCov)
}

#' @keywords internal
#' Prepare the matrices of focal cell type
prepare_focal <- function(sce,resolution,markerGene,focal_cell_index){
  ### Select cells of the focal cell type
  sce = sce[,focal_cell_index]
  resolution = min(resolution, ncol(sce))

  ### Divide cells along the cell trajectory into cell state bins
  sce$cell_state = Partition_cell_trajectory(sce$cell_trajectory,nbins = resolution)

  ### Compute the mean expression profile of cell state bin
  ExpBin = Seurat::AverageExpression(sce[markerGene,],group.by = 'cell_state',slot = 'counts',assays = 'RNA')[[1]]
  ExpBin = ExpBin[,order(readr::parse_number(colnames(ExpBin)))]

  ### Compute the median state-score(pseudo-time) cell state bin
  TimeBin = aggregate(sce$cell_trajectory,by=list(sce$cell_state),FUN=median)[,-1]
  TimeBin = sort(TimeBin,decreasing = FALSE)

  ### Extract the expression profile of individual cells
  cell_trajectory = as.matrix(sort(sce$cell_trajectory,decreasing = FALSE))
  ExpCell = GetAssayData(object = sce[markerGene,], slot = "counts",assay = 'RNA')[,rownames(cell_trajectory)]
  ExpCell = as.matrix(ExpCell)

  focal_input = list('TimeBin'=TimeBin,'ExpBin'=ExpBin,'ExpCell'=ExpCell,"cell_trajectory"=cell_trajectory)
  return(focal_input)
}

#' @keywords internal
#' Smooth the estimated cell-state abundance
smooth_abundance <- function(abundance,smoothMethod,neighbor,span,TimeBin){

  if(smoothMethod=='loess'){
  ### Use the loess method to smooth the estimated cell-state abundance
    abundance = apply(abundance,2,function(y){predict(loess(y~TimeBin,span = span ))})
  }
  else{
  ### Use the neighbor mean to smooth the estimated cell-state abundance
    num = min(neighbor,round(length(TimeBin)*0.2))
    neighbor_cell = sapply(1:length(TimeBin),function(i){order(abs(TimeBin[i]-TimeBin),decreasing = F)[1:num]})

    abundance = apply(abundance,2,function(y){
      sapply(1:ncol(neighbor_cell),function(i){mean(y[neighbor[,i]])})
    })
  }
  return(abundance)
}

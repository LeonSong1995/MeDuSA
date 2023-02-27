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


#' @keywords internal
#' Measure cel-state abundance from the scRNA-seq data 
Measure_abundance <- function(pseudotime,sampleID,resolution){
  pseudotime = pseudotime[names(sampleID)]
  bin = paste0('bin',cut(pseudotime,resolution,labels = FALSE, include.lowest = TRUE))
  breaks = aggregate(pseudotime,by=list(bin),FUN=min)
  breaks = breaks[order(readr::parse_number(breaks[,1])),2]
  breaks[1]=-Inf;breaks[length(breaks)+1]=Inf

  abundance_expect = sapply(unique(sampleID),function(id){
    pseudotime_temp = sort(pseudotime[which(sampleID == id)])
    pseudotime_temp = pseudotime_temp
    #count the cell number for each cell-state bin
    count_temp = sapply(2:length(breaks), function(currBreakIndex) {
      length(which(pseudotime_temp >= breaks[currBreakIndex-1] & pseudotime_temp < breaks[currBreakIndex]))
    })
    #normalize to the fractional abundance
    abundance_temp  =  count_temp/sum(count_temp)
    return(abundance_temp)
  })
  rownames(abundance_expect) = unique(bin)[order(readr::parse_number(unique(bin)))]

  return(abundance_expect)
}

#' @keywords internal
#' Evaluate the method performance
evalualtion <- function(abundance_simu,abundance_est){
  abundance_get = abundance_simu*0
  common_bin = intersect(rownames(abundance_get),rownames(abundance_est))
  common_ID = intersect(colnames(abundance_get),colnames(abundance_est))
  abundance_get[common_bin,common_ID] = abundance_est[common_bin,common_ID]

  #Normalize 
  abundance_get = sweep(abundance_get,2,colSums(abundance_get),'/')

  #1)-CCC
  CCC = sapply(colnames(abundance_simu),function(id){
    real = abundance_simu[,id]
    est = abundance_get[,id]
    c(DescTools::CCC(real,est)$rho.c[1])
  })
  CCC = unlist(CCC)
  #2)-Pearson R
  R = diag(cor(abundance_simu,abundance_get[,colnames(abundance_simu)]))
  #3)-RMSD
  RMSD = sqrt(colMeans((abundance_simu-abundance_get)^2))

  metric = rbind(CCC,R,RMSD)
  colnames(metric) = common_ID
  return(metric)
}

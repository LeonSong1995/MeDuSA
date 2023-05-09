##########################################################################################################################################################################################################################################################################################
##README: Functions for selecting marker genes
##########################################################################################################################################################################################################################################################################################
### MeduSA provides two methods to select marker genes over the given cell-state trajectory.
### ----------------------------------------------------------------------------------------------------------------------
### The first one is the Wilcoxon test: 
### MeduSA first divides the cells in the trajectory into a specified number of bins (by default, 10 bins are used). 
### For each bin, MeduSA applies the Wilcoxon-test implemented in the Seurat::FindMarkers function.
### The Wilcoxon test is used to compare gene expression levels between two groups of cells and determine whether the difference in expression is statistically significant. 
### In this case, the two groups of cells are the cells in the current bin being tested and all the other cells in the trajectory. 
### The Wilcoxon test is performed for each gene, and genes with significant differential expression are identified as marker genes for that particular bin.
### By applying the Wilcoxon test to each bin along the cell-state trajectory, MeduSA can identify marker genes that are specific to each stage of the trajectory.
### ----------------------------------------------------------------------------------------------------------------------
### The second one is the gam-wald test: 
### MeduSA associates genes along the cell-state trajectory using the generalized additive model (gam). 
### Only genes with a false discovery rate (FDR) adjusted p-value less than 0.01 are considered. 
### These significant genes are then ranked based on their association strength, which allows for the identification of the most relevant genes that are associated with the cell-state trajectory.
### To prevent certain cell states from being overrepresented, MeduSA divides the cell-state trajectory into a specified number of intervals (by default, 10 intervals are used). 
### Each gene is then assigned to the interval in which it has the highest mean expression. 
### For each interval, a set of top informative genes is selected as signature genes.
### ----------------------------------------------------------------------------------------------------------------------
### Once the gene selection process is complete, MeduSA then checks the representativeness of the selected genes for each cell-state bin. 
### This is done to ensure that the selected genes are indeed informative for the cell states they were selected for and are not biased towards a particular subset of cells. 
### ----------------------------------------------------------------------------------------------------------------------

#' @keywords internal
#' Select marker genes to cover the cell trajectory
MeDuSA_marker <- function(sce,bulk,geneNumber,nbins,family,k,ncpu,method){

  #1)---divides cells to cell-state bins
  cellStateBin = Partition_cell_trajectory(sce$cell_trajectory, nbins)
  Idents(sce) = cellStateBin
  
  #2)---Select common gene expressed in both scRNA-seq and bulk RNA-seq data
  commonGene = rownames(sce)[rownames(sce) %in% rownames(bulk)]

  #2.a)---Select marker genes (method-1: using Wilcox test)
  if(method=='wilcox'){
    mk_list = MK_wilcox(sce[commonGene,],cellStateBin,ncpu=ncpu,geneNumber=geneNumber)
  }

  #2.b)---Select marker genes (method-2: using Gam-Wald test)
  if (method=='gam'){
    mk_list = MK_gam(sce[commonGene,],cellStateBin,ncpu=ncpu,family=family,k=k, geneNumber=geneNumber)
  }

  #3)---Along the cell trajectory, check the coverage of marker genes
  markerGene = MK_balance(sce,bulk,mk_list)

  return(markerGene)
}

#' @keywords internal
#' Along the cell trajectory, partition cells into cell state bins.
Partition_cell_trajectory <- function(XY,nbins){
  cellLocation = cut(XY, nbins, labels = FALSE, include.lowest = TRUE)
  cellLocation = paste0('bin',cellLocation)
  return(cellLocation)
}

#' @keywords internal
#' Select marker genes using wilcox
MK_wilcox <- function(sce,cellStateBin,ncpu,geneNumber){

  mt_gene = grep(pattern = "MT-",rownames(sce),value = F)
  rp_gene = grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",rownames(sce),value = F)
  sce = sce[-c(mt_gene,rp_gene),]

  ### To get a robust result, we agrregate the cell-state bins with cell number <= 20
  smallBin = names(which(table(cellStateBin) <= 20))
  stateBin = aggregate(sce$cell_trajectory,by=list(cellStateBin),FUN=median)
  binName = stateBin[,1];stateBin = stateBin[,-1]
  names(stateBin) = binName
  for(bin in smallBin){
    mergeBin = names(which.min(abs(stateBin[bin] - stateBin[!names(stateBin) %in% smallBin])))
    cellStateBin[which(cellStateBin==bin)] = mergeBin
  }
  Idents(sce) = cellStateBin
  # print(table(Idents(sce)))

  ### Register CPU cores
  ncpu = min(ncpu,parallel::detectCores())
  message('\n',paste0(paste0('Select genes using Wilcoxon test with ',ncpu)),' cores.')
  cl = parallel::makeCluster(ncpu)

  ### Register environment
  GeneNumberBin = round(geneNumber/length(unique(cellStateBin)))
  parallel::clusterExport(cl=cl, varlist=c("sce","cellStateBin","GeneNumberBin"), envir=environment())
  doSNOW::registerDoSNOW(cl)

  ### Show the process
  binName = names(which(table(cellStateBin)>3))
  binNum = length(binName)
  pb = utils::txtProgressBar(min = 1, max = binNum, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`

  ### Select genes (Wilcox Test)
  mk = foreach::foreach(BinId = 1:binNum, .options.snow = opts) %dopar2% {
    mk_focal = Seurat::FindMarkers(sce,ident.1 = binName[BinId],verbose = F,only.pos = T,
                                   slot = "data",assay = sce@active.assay,logfc.threshold = 0.3)
    mk_focal = mk_focal[mk_focal$p_val_adj < 0.01,]

    ## Order the gene with p_val=0 using fold changes
    index_p0 = which(mk_focal$p_val == 0)
    index_p0 = index_p0[order(mk_focal[index_p0, 'avg_log2FC'], decreasing = TRUE)]
    index = setdiff(seq_len(nrow(mk_focal)), index_p0)
    index = c(index_p0, index)
    mk_focal = mk_focal[index,]
    return(mk_focal)
  }
  parallel::stopCluster(cl)
  names(mk) = binName

  ### Compute the mean expression of cell-state bin
  expMean = Seurat::AverageExpression(sce,slot = 'data',assays = sce@active.assay)[[1]]
  mk_filter = sapply(binName,function(bin){
    gene = mk[[bin]]
    gene = gene[which(colnames(expMean)[apply(expMean[rownames(gene),],1,which.max)]==bin),]
    return(list(gene))
  })

  ### Select top genes
  if(min(sapply(mk_filter,nrow)) <= 5){
    warning(paste0(paste0('The number of signature genes in ',
              paste(names(which(sapply(mk_filter,nrow)<5)),collapse =', ')),' is smaller than 5!'))
    mk = sapply(mk_filter,function(gene){list(rownames(gene)[1:GeneNumberBin])})
  }else{
    mk = sapply(mk,function(gene){list(rownames(gene)[1:GeneNumberBin])})
  }

  return(mk)
}


#' @keywords internal
#' Select marker genes using gam
MK_gam <- function(sce,cellStateBin,ncpu,family,k,geneNumber){

  mt_gene = grep(pattern = "MT-",rownames(sce),value = F)
  rp_gene = grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",rownames(sce),value = F)
  sce = sce[-c(mt_gene,rp_gene),]

  ### Register CPU cores
  ncpu = min(ncpu,parallel::detectCores())
  message('\n',paste0(paste0('Select genes using Gam-Wald test with ',ncpu)),' cores.')
  cl = parallel::makeCluster(ncpu)

  ### Register environment
  eligibleGene = names(which((Matrix::rowSums(sign(sce@assays$RNA@counts)) / ncol(sce)) > 0.1))
  sce = sce[eligibleGene,]
  exprsData = GetAssayData(object = sce, slot = "data",assay = sce@active.assay)
  space = sce$cell_trajectory
  parallel::clusterExport(cl=cl, varlist=c("exprsData","space","k"),envir=environment())
  doSNOW::registerDoSNOW(cl)

  ### Show the process
  pb = utils::txtProgressBar(min = 1, max = nrow(exprsData), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`

  ### Compute gene association strength (gam-wald Test)
  geneId = NULL
  Chi = foreach::foreach(geneId = 1:nrow(exprsData), .options.snow = opts) %dopar2% {
    gam_mod = mgcv::gam(exprsData[geneId,] ~ 1+s(space,k=k,bs='cr',fx=FALSE),
                        family = family,method='GCV.Cp')
    mgcv::anova.gam(gam_mod)$chi.sq
  }
  parallel::stopCluster(cl)
  Chi  = unlist(Chi)

  ### Compute FDR p-values
  P_adj = p.adjust(pchisq(Chi,df=1,lower.tail = F))
  names(P_adj) = names(Chi) = rownames(sce)

  ### Compute the mean expression of cell-state bin
  expMean = Seurat::AverageExpression(sce,slot = 'counts',assays = 'RNA')[[1]]
  maxBin = colnames(expMean)[apply(expMean,1,which.max)]
  maxValue = apply(expMean,1,max)
  Info = data.frame('bin' = maxBin,'chi' = Chi,'p_val_adj' = P_adj,'avg_exp'=maxValue)
  Info = Info[Info$p_val_adj < 0.01,]

  ### Select genes for each cell state bin
  GeneNumberBin = round(geneNumber/length(unique(cellStateBin)))
  mk = sapply(unique(cellStateBin),function(focalBin){
    Info_focal = Info[Info$bin==focalBin,]
    Info_focal = Info_focal[order(Info_focal$chi,decreasing = T),]
    mk_focal = rownames(Info_focal)[1:GeneNumberBin]
    list(mk_focal)
  })

  return(mk)
}

#' @keywords internal
#' Check gene coverage
MK_balance <- function(sce,bulk,mk_list){

  ### Compute the mean maker expression of each cell state bin
  mk = Reduce(intersect,list(unique(unlist(mk_list)),rownames(sce),rownames(bulk)))
  expMean = Seurat::AverageExpression(sce[mk,],slot = 'counts',assays = "RNA")[[1]]

  ### Re-select genes (to ensure the balance coverage)
  # geneNumberBin = round(min(table(apply(expMean,1,which.max)))*1)
  geneNumberBin = round(mean(table(apply(expMean,1,which.max))))
  MarkerGene = lapply(mk_list,function(g){na.omit(g[g %in% mk][1:geneNumberBin])})
  MarkerGene = unique(unlist(MarkerGene))

  ### Remove the over-weighted genes (> mean+ 3*sd)
  expMean_marker = apply(expMean[MarkerGene,],1,mean)
  MarkerGene = names(which(expMean_marker < (mean(expMean_marker) + 3*sd(expMean_marker))))

  return(MarkerGene)
}

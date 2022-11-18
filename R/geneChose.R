#' @keywords internal
#Group cells into cell-state-bins
cluster =  function(XY,nbins){
  left = min(XY)-1e-7
  right = max(XY)+1e-7
  breaks = seq(left,right,(right-left)/nbins)
  cellLocation = rep(NA,length(XY))
  for(i in 1:length(breaks)-1){cellLocation[XY >= breaks[i] & XY < breaks[i+1]]=i}
  return(cellLocation)
}

#' @keywords internal
#Select genes via GAM
geneAsso = function(space,exprsData,CellBin,maxgene,cov=NULL,family,k,ncpu){

  #1)---Pre-adjust the covariates
  if(!is.null(cov)){
    message('Adjust covariates of gene expression')
    formula = as.formula(paste0('y~',paste(colnames(cov),collapse='+')))
    exprsData = t(apply(exprsData,1,function(y){
      dat = data.frame(y=y,cov)
      residuals(lm(formula=formula,dat))
    }))
  }

  #2)---Cell-state-bin expression profile
  groupRef = t(aggregate(t(exprsData),by=list(CellBin),FUN=mean)[,-1])
  colnames(groupRef) = paste0(unique(CellBin),'c')
  genes_to_take = rownames(groupRef)
  cellNum = ncol(exprsData)

  #3)---RUN GAM to associate gene expression with the cell trajectory
  ncpu = min(ncpu,detectCores())
  cl = parallel::makeCluster(ncpu)

  parallel::clusterExport(cl=cl, varlist=c("exprsData","space","k"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(min = 1, max = nrow(exprsData), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  geneNumber = NULL
  fitF= foreach::foreach(geneNumber = 1:nrow(exprsData), .options.snow = opts) %dopar2% {
    gam_mod=mgcv::gam(exprsData[geneNumber,] ~ 1+s(space,k=k,bs='cr',fx=FALSE),family = family,method='GCV.Cp')
    gam_mod=mgcv::anova.gam(gam_mod)$chi.sq
    gam_mod
  }
  parallel::stopCluster(cl)

  fitF  = unlist(fitF)
  names(fitF) = genes_to_take
  res = data.frame(colnames(groupRef)[apply(groupRef,1,which.max)],fitF)
  colnames(res)= c("group", "Chi_score")
  rownames(res)= rownames(exprsData)

  #4)---Rank genes based on their Chi_score
  res = res[order(res$Chi_score,decreasing = T),]

  #5)---Select the top N genes
  geneNumberCluster = ceiling(maxgene/ncol(groupRef))
  topGenes = unlist(sapply(colnames(groupRef),function(bin){
    index = which(as.character(res$group)==as.character(bin))
    temp = rownames(res[index,])
    na.omit(temp[1:min(geneNumberCluster,length(temp))])
  }))

  #6)---Remove the over-weighted genes (>mean+3sd)
  expM = apply(exprsData[topGenes,],1,mean)
  topGenes=topGenes[-which(expM>(mean(expM)+2*sd(expM)))]

  ###End
  return(list('bestGenes'=topGenes,'Chi'=fitF))
}


#' @keywords internal
geneSelect = function(exprsData,space,bulk,maxgene,nbins=10,cov=NULL,family='gaussian',k=10,ncpu){

  message('\n',paste0(paste0('Select genes with ',ncpu)),' cores.')
  #1)---Choose genes both expressed in scRNA-seq and bulk RNA-seq
  commGene = intersect(rownames(exprsData),rownames(bulk))
  exprsData = exprsData[commGene,]

  #2)---Choose genes expressed in >10% cells
  eligible_gene = names(which(rowSums(sign(exprsData))/ncol(exprsData)>0.1))
  exprsData = exprsData[eligible_gene,]
  print(dim(exprsData))

  #3)---Assign cells to cell-trajectory bins
  CellBin = cluster(space, nbins)

  #4)---Choose genes using the GAM
  gw = geneAsso(space=space,exprsData=exprsData,CellBin=CellBin,maxgene=maxgene,cov=cov,family=family,k=k,ncpu=ncpu)
  g = gw$bestGenes
  chi = gw$Chi

  return(list('g'=g,'chi'=chi))
}

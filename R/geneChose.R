
#' @keywords internal
#Group cells into cell-bins
cluster =  function(XY,nbins){
  left = min(XY)-1e-7
  right = max(XY)+1e-7
  breaks = seq(left,right,(right-left)/nbins)
  cellLocation = rep(NA,length(XY))
  for(i in 1:length(breaks)-1){
    cellLocation[XY >= breaks[i] & XY < breaks[i+1]]=i
  }
  return(cellLocation)
}

#' @keywords internal
VariableGenes = function(gene,ratio){
  value = length(which(gene > min(gene)))
  cellnum = length(gene)
  if((value/cellnum)>ratio){TRUE}else{FALSE}
}

#' @keywords internal
#Select genes via GAM
geneAsso = function(space,exprsData,TJB,maxgene,cov=NULL,family,k,ncpu){  
  
  #TJB-level profile
  groupRef = t(aggregate(t(exprsData),by=list(TJB),FUN=mean)[,-1])
  colnames(groupRef) = paste0(unique(TJB),'c')
  genes_to_take = rownames(groupRef)
  cellNum = ncol(exprsData)
  
  #Prepare the covariates
  U = rep(1,cellNum)
  if(!is.null(cov)){U = cov}
  
  #RUN GAM
  ncpu = min(ncpu,detectCores())
  cl = parallel::makeCluster(ncpu)
  parallel::clusterExport(cl=cl, varlist=c("U","exprsData","space","k"),
                          envir=environment())
  doSNOW::registerDoSNOW(cl)
  pb = utils::txtProgressBar(min = 1, max = nrow(exprsData), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  `%dopar2%` = foreach::`%dopar%`
  geneNumber = NULL
  fitF= foreach::foreach(geneNumber = 1:nrow(exprsData), .options.snow = opts) %dopar2% {
    gam_mod=mgcv::gam(exprsData[geneNumber,] ~ U+s(space,k=k,bs='cr',fx=FALSE),family = family,method='GCV.Cp')
    gam_mod=mgcv::anova.gam(gam_mod)$chi.sq
    gam_mod
  }
  parallel::stopCluster(cl)

  # # *Run with polynomial regression###
  # ncpu = min(ncpu,detectCores())
  # cl = parallel::makeCluster(ncpu)
  # X=poly(space,k)
  # parallel::clusterExport(cl=cl, varlist=c("X","exprsData","U"),
  #                         envir=environment())
  # doSNOW::registerDoSNOW(cl)
  # pb = utils::txtProgressBar(min = 1, max = nrow(exprsData), style = 3)
  # progress = function(n) setTxtProgressBar(pb, n)
  # opts = list(progress = progress)
  # `%dopar2%` = foreach::`%dopar%`
  # geneNumber = NULL
  # fitF = foreach::foreach(geneNumber = 1:nrow(exprsData), .options.snow = opts) %dopar2% {
  #   unlist(summary(aov(lm(exprsData[geneNumber,]~X+U))))['F value1']
  # }
  # parallel::stopCluster(cl)

  fitF  = unlist(fitF)
  names(fitF) = genes_to_take
  res = data.frame(colnames(groupRef)[apply(groupRef,1,which.max)],fitF)
  colnames(res)= c("group", "Chi_score")
  rownames(res)= rownames(exprsData)
  
  #Rank genes based on their Chi_score
  res = res[order(res$Chi_score,decreasing = T),]
  
  #Select the top N genes
  geneNumberCluster = ceiling(maxgene/ncol(groupRef))
  topGenes = unlist(sapply(colnames(groupRef),function(bin){
    index = which(as.character(res$group)==as.character(bin))
    temp = rownames(res[index,])
    na.omit(temp[1:min(geneNumberCluster,length(temp))])
  }))
  
  #Remove the over-weighted genes (>mean+2sd)
  expM = apply(exprsData[topGenes,],1,mean)
  topGenes=topGenes[-which(expM>(mean(expM)+2*sd(expM)))]
  
  
  return(list('bestGenes'=topGenes,'Chi'=fitF))
}


#' @keywords internal
geneSelect = function(exprsData,space,bulk,maxgene,nbins=10,cov=NULL,family='gaussian',k=10,ncpu){
  
 message('\n',paste0(paste0('Select genes with ',ncpu)),' cores.')
  #Choose genes both expressed in scRNA-seq and bulk RNA-seq
  commGene = intersect(rownames(exprsData),rownames(bulk))
  exprsData = exprsData[commGene,]
  
  #Choose cell-expressed genes
  exprsData = exprsData[apply(exprsData,1,VariableGenes,0.1),]
  print(dim(exprsData))
  
  #TJB = TraJectory Bins (assign cells to cell-trajectory bins) 
  TJB = cluster(space, nbins)
  
  #Choose genes with the generalized nonlinear additive model 
  gw = geneAsso(space=space,exprsData=exprsData,TJB=TJB,maxgene=maxgene,cov=cov,family,k=k,ncpu=ncpu)
  g = gw$bestGenes
  chi = gw$Chi

  return(list('g'=g,'chi'=chi))
}


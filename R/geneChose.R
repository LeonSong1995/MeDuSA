
#############################genes for cell-type level deconvolution ############################################################################
#' @keywords internal
SignatureGenerator = function(sce){

	cluster.markers = matrix(0,1,6)
	cluster.markers = as.data.frame(cluster.markers)
	colnames(cluster.markers) = c("p_val", "avg_logFC","pct.1", "pct.2", "p_val_adj", "cluster" )
	cluster.markers = cluster.markers[-c(1),]

	Idents(sce) = sce$cellType

	for (i in unique(sce$cellType)){
		print(paste('Finding signature genes for:',i,sep=""))
		cluster1.markers = FindMarkers(sce, ident.1 = i, min.pct = 0.25)
		temp = as.data.frame(matrix(0,nrow(cluster1.markers),1))
		cluster1.markers$cluster = temp[,1]
		cluster1.markers = cluster1.markers[which(cluster1.markers$avg_logFC>0 & cluster1.markers$p_val<0.01),]
		#Select the top 100 genes for each cluster
		if(nrow(cluster1.markers)>100){cluster1.markers = cluster1.markers[1:100,]}
		cluster.markers = rbind(cluster.markers,cluster1.markers)
	}

	FinalMarkers<-rownames(cluster.markers)
	return (FinalMarkers)

}



#############################genes for single-cell level deconvolution ############################################################################
#' @keywords internal
Initial = function(XY, refNames,nbin){
	#Add grid
	nbins = round(length(refNames)*nbin,0)
	currXY = XY
	breaks = seq(min(currXY)-10^-7,max(currXY)+10^-7, (max(currXY)-min(currXY)+2*10^-7)/ceiling(nbins))
	grid = rep(NA,ceiling(nbins))
	cellLocationOnGrid = rep(NA,length(currXY))
	for(currBreakIndex in 2:length(breaks)){
		cellLocationOnGrid[which(currXY>breaks[currBreakIndex-1] & currXY<breaks[currBreakIndex])] = currBreakIndex-1
	}
	tab = table(cellLocationOnGrid)
	
	#relocation (Merge small grid)
	# for (i in as.numeric(names(tab)[which(tab<2)])){
 #    if (i==min(cellLocationOnGrid)){
 #      cellLocationOnGrid[which(cellLocationOnGrid==i)] = as.numeric(names(tab))[which(i==as.numeric(names(tab)))+1]
 #    }else{
 #      cellLocationOnGrid[which(cellLocationOnGrid==i)] = as.numeric(names(tab))[which(i==as.numeric(names(tab)))-1] 
 #    }
 #    tab = table(cellLocationOnGrid)
 #  }

	grid[as.numeric(names(tab))] = tab
	gridCellsToUse = which(!is.na(grid))

	#Choose grid center
	chosenCellsForCluster = unlist(lapply(gridCellsToUse, function(currCell){
	chosenCell = which(cellLocationOnGrid==currCell)
	if(length(chosenCell)>1){
		center = log2(mean(XY[chosenCell,])+1)
		chosenCell = chosenCell[which.min(fields::rdist(t(as.matrix(center)),log2(currXY[chosenCell,]+1)))]
	}
	chosenCell
	}))
	chosenCellsForCluster
}

#' @keywords internal
geneANOVA = function(XY,exprsData,centerSelect,refNames,op,maxgene){	
	#K-Means
	mode = paste(unique(refNames),kmeans(XY,centers = XY[centerSelect,])$cluster,sep='-')
	names(mode) = rownames(XY)
	center = as.matrix(kmeans(XY,centers = XY[centerSelect,])$centers)
	rownames(center) =  names(table(mode))

	#Re-group
	Regroup = unlist(sapply(1:length(table(mode)),function(id){
		modeChose = names(table(mode))[id]
		currSpace = as.matrix(XY[which(mode==modeChose),])
		currName = rownames(currSpace)
		Num = ceiling(nrow(currSpace)*op)
		currCenter = as.matrix(center[id,])
		cellDist = as.vector(fields::rdist(t(as.matrix(currCenter)),currSpace))
		names(cellDist) = currName
		out = names(sort(cellDist,decreasing = T)[1:Num])
		regroup = sapply(1:length(out),function(j){
			remap = as.vector(fields::rdist(t(as.matrix(currSpace[out[j],])),center[-id,]))
			names(remap) = names(table(mode))[-id]
			names(sort(remap,decreasing = F)[1])
		})
	  names(regroup) = out
	  regroup
	}))
	modeBig = c(mode,Regroup)
	chosenNeigList = sapply(names(table(modeBig)),function(id){names(which(modeBig ==id))})

	#Cell-level profile
	cellRef = do.call(cbind,sapply(seq(1:length(chosenNeigList)),function(cell){
		exprsData[,chosenNeigList[[cell]]]
	}))
	colnames(cellRef) = unlist(sapply(seq(1:length(chosenNeigList)),function(cell){
		paste(cell,rep('c',length(chosenNeigList[[cell]])),sep='')
	}))
	cellNamesForAnova = as.vector(colnames(cellRef))

	#Cluster-level profile
	groupRef = sapply(seq(1,length(chosenNeigList)),function(cell){
		temp = exprsData[,chosenNeigList[[cell]]]
		if(!is.null(nrow(temp))){
			return(rowMeans(temp))
		}else{return(temp)}
	})
	colnames(groupRef) = names(table(mode))
	
	#ANOVA
	genes_to_take = rownames(groupRef)
	dat = cbind(rep(0,length(genes_to_take)),cellRef[genes_to_take,])
	group = c("",cellNamesForAnova)
	dmat <- stats::model.matrix(~ group)
	fit <- limma::eBayes(limma::lmFit(dat, dmat)[,-1])
	fitF = fit$F
	res = as.data.frame(cbind(gsub("group","",colnames(fit$coefficients)[ apply(fit$coefficients,1,which.max)]),fitF))
	colnames(res) = c("group", "score")
	allGenes = apply(as.data.frame(unique(cellNamesForAnova)),1,function(bin){
		Indexes = which(as.character(res$group)==as.character(bin))
		(genes_to_take[Indexes])[order(res$score[Indexes],decreasing = T)]
	})

	#Choose gene 
	bestKappa = Inf
	bestG = 0
	mul = 1
	maxNumberOfGenesPerCell = round(maxgene/ncol(groupRef))
	minNumberOfGenesPerCell = 2
	bestGenes = c()
	indexRange = minNumberOfGenesPerCell:maxNumberOfGenesPerCell
	for (i in indexRange){
	  selectedGenes = unique(as.character(unlist(lapply(1:length(allGenes), function(listIndex){
		unlist(allGenes[listIndex])[which(!is.na(unlist(allGenes[listIndex])[1:as.numeric(i*mul)]))]
	  }))))
	  currcellRef = exprsData[match(selectedGenes, row.names(exprsData)),]
	  currgroupRef = groupRef[match(selectedGenes, row.names(groupRef)),]
	  #ecor: a covariance function written in C++ (see REML part)
	  newKappa = kappa(ecor(currgroupRef),exact = T)
	  
	  if (newKappa < bestKappa){
		bestKappa = newKappa
		bestG = i
		bestGenes = unique(selectedGenes)
	  }
	}
	
	
	return(list('bestKappa'=bestKappa,'bestGenes'=bestGenes))
}

#' @keywords internal
checkVariableGenes = function(a, ratio) {
	count_nonZeros = length(which(a > min(a)))
	if (count_nonZeros/length(a) > ratio) {var(a)/ mean(a)} else {0}
}

#' @keywords internal
geneSelect = function(exprsData,cellType,space,bulk,op,maxgene,nbin){
	commgene = intersect(rownames(exprsData),rownames(bulk))
	#Choose cell-expressed genes
	geneVarianceRef = apply(exprsData,1,function(gene){checkVariableGenes(as.numeric(as.matrix(gene)),0.1)})
  geneVarianceFinalRef = sort(geneVarianceRef[geneVarianceRef>0],decreasing = T)
	mutualGenes = intersect(names(geneVarianceFinalRef),commgene)
	exprsData = exprsData[mutualGenes,]

	#Run ANOVA
	InitialCenter = Initial(space,cellType,nbin)
	gw = geneANOVA(space,exprsData,InitialCenter,cellType,op,maxgene)
	g = gw$bestGenes
	
	return(g)
}

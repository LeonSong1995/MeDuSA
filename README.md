# MCTD: Mixed model-based Cell Trajectory Deconvolution. 
**Author: Liyang Song <songliyang@westlake.edu.cn>**    


## Description
MCTD is a fine-resolution deconvolution method used to predict cell abundance along a predefined cell trajectory in the given bulk RNA-seq data.


## Installation
```R
install.packages("devtools")

##Please install the "Seurat" first. (https://satijalab.org/seurat/)
install.packages("Seurat")
library(Seurat)

##R version need > 3.5.0
devtools::install_github("LeonSong1995/MCTD", build_vignettes=F)
```


## Usage
The function for **Cell-Trajectory** deconvolution in this package is `MCTD`. It needs:  
1. Bulk RNA-seq data.  A matrix of bulk RNA-seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
2. Single-cell RNA-seq data. A [Seurat](https://satijalab.org/seurat/) obejct of the scRNA-seq data. 

## Example
MCTD provides test data to show how to use.
```R
library(MCTD)
##When you load MCTD, you will get these attached test data:
data(ref)
data(cellType)
data(cellTrjaectory)
data(bulk)

##You need to format them into the 'Seurat' obejct:
sce = CreateSeuratObject(ref)
sce$cellType = cellType
sce$cellTrjaectory = rep(0,ncol(sce))
sce$cellTrjaectory[rownames(Trajectory)]=Trajectory

## Now, you can run cell-trajectory deconvolution.
##Run MCTD (with 6 CPU cores)
CellAbundance = MCTD(bulk=bulk,sce=sce,selectCellType='Epithelium',ncpu=6)$abundance

##Detailed tutorials for the MCTD can be found via: 
help(MCTD)
```

## Errors
The following errors mean that the mixed model can not converge on your data. When you see them, you may need to check the signature genes or the reference scRNA-seq data.  
`REML ERROR!:V matrix is not positive.`  
`REML ERROR!: the X^t * V^-1 * X matrix is not invertible,please check the Signature Genes.`


## Contact
If you have any questions for MCTD, please create an issue here or contact us <songliyang@westlake.edu.cn>.



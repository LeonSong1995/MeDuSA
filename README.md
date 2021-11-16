# MeDuSA: mixed model-based deconvolution of cell-state abundance.

**Author: Liyang Song <songliyang@westlake.edu.cn>**    


## Description
MeDuSA is a fine-resolution cellular deconvolution method, with the aim to use reference scRNA-seq data to predict cell abundance distributed along a cell-state trajectory in a bulk RNA-seq data. MeDuSA is well-suited for biological scenarios in which the underlying mechanisms are associated with continuous transitions of cell-states.


## Installation
```R
install.packages("devtools")

##Please install the "Seurat" first. (https://satijalab.org/seurat/)
install.packages("Seurat")
library(Seurat)

##R version need > 3.5.0
devtools::install_github("LeonSong1995/MeDuSA", build_vignettes=F)
```


## Usage
The function for **cell-state abundance** deconvolution in this package is `MeDuSA`, which needs:  
1. Bulk RNA-seq data.  A matrix of bulk RNA-seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
2. Single-cell RNA-seq data. A [seurat](https://satijalab.org/seurat/) obejct of the reference scRNA-seq data. 

## Example
MeDuSA package provides test data to show how to use.
```R
##Library the package
library(MeDuSA)

##Load the attached test data
data(ref)
data(cellType)
data(Trajectory)
data(bulk)

##Build the seurat obejct
sce = CreateSeuratObject(ref)
sce$cellType = cellType
sce$cellTrajectory = rep(0,ncol(sce))
sce$cellTrajectory[rownames(Trajectory)]=Trajectory

##Run MeDuSA (run with 6 courses):
CellStateAbundance = MeDuSA(bulk=bulk,sce=sce,selectCellType='Epithelium',ncpu=6)
Abundance = CellStateAbundance$abundance
SignatureGene = CellStateAbundance$gene
PseudoTime = CellStateAbundance$PesudoTimeCellbin

##Documents
help(MeDuSA)
```

## Errors
The following errors mean that the MeDuSA model can not converge on your data. When you see them, you may need to check the signature genes or the reference scRNA-seq data.  
`REML ERROR!:V matrix is not positive.`  
`REML ERROR!: the X^t * V^-1 * X matrix is not invertible,please check the Signature Genes.`


## Contact
If you have any questions for MeDuSA, please create an issue here or contact the author <songliyang@westlake.edu.cn>.   

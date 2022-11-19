# MeDuSA: mixed model-based deconvolution of cell-state abundance.

**Author: Liyang Song <songliyang@westlake.edu.cn>**    


## Description
![Image text](https://github.com/LeonSong1995/MeDuSA/blob/master/schematic/schematic.jpg)

**MeDuSA** is a fine-resolution cellular deconvolution method that leverages scRNA-seq data as a reference to estimate cell-state abundance in bulk RNA-seq data. **MeDuSA** features the use of a linear mixed model (LMM) to fit a cell state in question (either a single cell or the mean of multiple cells) as a fixed effect and the remaining cells of the same cell type individually as random effects accounting for correlations between cells. This model improves the deconvolution accuracy because the random-effect component allows each cell has a specific weight on bulk gene expression, resulting in a better capturing of variance in bulk gene expression. This model aslo ameliorates the collinearity problem between cells at the focal state (fitted as a fixed effect) and those at adjacent states (fitted as random effects) because of the shrinkage of random effects.

## Installation
```R
#1)---Please install the Seurat. (https://satijalab.org/seurat/)
install.packages("Seurat")

#2)---Please install the BiocParallel. (https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BiocParallel")

#3)---Install the MeDuSA (R version > 3.5.0)
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("LeonSong1995/MeDuSA", build_vignettes=F)
```


## Usage
The function for **cell-state abundance** deconvolution in this package is `MeDuSA`, which needs:  
1. Bulk RNA-seq data.  A matrix of bulk RNA-seq data. Each row corresponds to a specific gene and each column corresponds to a particular sample.
2. Single-cell RNA-seq data. A [Seurat](https://satijalab.org/seurat/) obejct of the reference scRNA-seq data. 

## Example
MeDuSA package provides test data to show how to use.
```R
#1)---Library the package
library(MeDuSA)

#2)---Load the test data
data(ref)
data(cellType)
data(Trajectory)
data(bulk)

#3)---Build the scRNA-seq reference (Seurat obejct)
sce = CreateSeuratObject(ref)
sce$cell_type = cellType
sce$cell_trajectory = rep(0,ncol(sce))
sce$cell_trajectory[rownames(Trajectory)]=Trajectory

#4)---Run MeDuSA (2 cores)
csab = MeDuSA(bulk=bulk,sce=sce,select.ct='Epithelium',ncpu=2)

#5)---Documents
help(MeDuSA)
```

## Errors
The following errors mean that the MeDuSA model can not converge on your data. You may need to check the signature genes or the reference scRNA-seq data when you see them.

`REML ERROR!:V matrix is not positive.`  
`REML ERROR!: the X^t * V^-1 * X matrix is not invertible,please check the Signature Genes.`


## Contact
If you have any questions for MeDuSA, please contact the author <songliyang@westlake.edu.cn>.   

---
layout: page
title: Example Analysis
description: ~
---
This tutorial offers an illustrative analysis of the human monocytes data from [Oetjen et al., 2018](https://insight.jci.org/articles/view/124928) using MeDuSA. Prior to running the analysis, it is important to ensure that the MeDuSA package has been installed. For installation instructions, please refer to the following [link](https://github.com/LeonSong1995/MeDuSA).


## Input Data
`MeDuSA` requires two types of input data:
- Bulk RNA-seq data. 
- Single cell RNA-seq (scRNA-seq) data, which should be provided in the form of a Seurat object that includes the annotated cell-state trajectory and cell types. 

For how to prepare the cell-state trajectory data, please read the section of `Prepare reference scRNA-seq data` in this tutorial. 

The input data required for running this tutorial can be downloaded from the following [link](https://github.com/LeonSong1995/MeDuSA). 
Detailed information regarding the input data is provided as follows.
### 1. Bulk RNA-seq data
```r
#### load the example bulk RNA-seq data, 
bulk = readRDS("./Monocytes_bulk.rds")
class(bulk)
"matrix" "array" 
```
The bulk RNA-seq data is represented in a matrix format, where each row corresponds to a specific gene and each column corresponds to a particular sample.

### 2. Reference scRNA-seq data
```r
#### load the example scRNA-seq data, 
sce = readRDS("./Monocytes_sce.rds")
class(sce)
[1] "Seurat"
attr(,"package")
[1] "SeuratObject"

sce$cell_trajectory[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                 0.9220488                  0.5167408                  0.4567616 
		 
sce$cell_type[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                     "mon"                      "mon"                      "mon" 
```
The reference scRNA-seq data must be in the Seurat object format, where the cell-state trajectory is stored in `sce$cell_trajectory` and the cell-type is stored in `sce$cell_type`. For further information about Seurat, please visit the following [link](https://satijalab.org/seurat/).



## Cell State Deconvolution
```r
library(MeDuSA)
help(MeDuSA)
``` 
### 1. Basic usage of MeDuSA
This section provides an introduction to the basic usage of MeDuSA.
- bulk: A matrix of bulk RNA-seq data. 
- sce: A Seurat object of scRNA-seq data.  
- select.ct: A character variable indicating the focal cell type.
- markerGene: A character vector containing the marker genes across the cell-state trajectory.If not provided, MeDuSA will utilize the `MeDuSA_marker` function to select marker genes for the analysis.
- resolution: A numeric variable used to specify the number of cell-state bins along the cell trajectory.
- smooth: A boolean variable to determine whether to smooth the estimated cell-state abundance.
- span: A numeric variable to control the degree of smoothing.
- fractional: A boolean variable to determine whether to normalize the estimated cell-state abundance to the fractional abundance (0-1).
- ncpu: The number of CPU cores to be used. 

For further details about the parameters, please refer to this [link](https://github.com/LeonSong1995/MeDuSA).
```r
MeDuSA_obj = MeDuSA(bulk,sce,
                  select.ct = 'mon',markerGene = NULL,span = 0.35,
		  resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)		 
```
The results are stored in MeDuSA_obj@Estimation
```r
#The estimated cell-state abundance
MeDuSA_obj@Estimation$cell_state_abundance[1:3,1:3]
               A           B           C
bin1 0.018230362 0.012866330 0.015188014
bin2 0.011331655 0.009618875 0.010931174
bin3 0.007009345 0.007740727 0.008308389

#The median state (pseudo-time) of cell-state bins
MeDuSA_obj@Estimation$TimeBin[1:3]
[1] 0.006934536 0.031808065 0.050630336

#The used marker genes
MeDuSA_obj@Estimation$markerGene[1:3]
[1] "FCGR3A" "IFITM2" "IFITM3"
```

### 2. How to select marker genes
MeDuSA allows users to input their own cultivated marker genes. Additionally, MeDuSA provides two methods for selecting marker genes that are representative of the cell-state trajectory.
- wilcox test: 
- gam-wald test:


```r
library(ggplot2)
abundance = MeDuSA_obj@Estimation$cell_state_abundance

```
### 3. How to include other cell types as covariates

### 4. How to use the mode of Conditional Autoregressive (CAR)

### 5. Do I need to normalize the data 


## Prepare reference scRNA-seq data

## Compare the estimated cell-state abundance to the expected truth



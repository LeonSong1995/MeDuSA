---
layout: page
title: Example Analysis (monocytes)
description: ~
---
This tutorial offers an illustrative analysis of the human monocytes data from [Oetjen et al., 2018](https://insight.jci.org/articles/view/124928) using MeDuSA. Prior to running the analysis, it is important to ensure that the MeDuSA package has been installed. For installation instructions, please refer to the following [link](https://github.com/LeonSong1995/MeDuSA).


## Input data
`MeDuSA` requires two types of input data:
- Bulk RNA-seq data. 
- Single cell RNA-seq (scRNA-seq) data, which should be provided in the form of a Seurat object that includes the annotated cell-state trajectory and cell types. 

For how to prepare the cell-state trajectory data, please read the section of `Prepare reference data` in this tutorial. 

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



## Cell-state deconvolution analysis
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
- wilcox test: MeduSA first divides the cells in the trajectory into a specified number of bins. For each bin, MeduSA applies the wilcox test implemented in the `Seurat::FindMarkers` function. The wilcox test is used to compare gene expression levels between two groups of cells and determine whether the difference in expression is statistically significant. In this case, the two groups of cells are the cells in the current bin being tested and all the other cells in the trajectory. The wilcox test is performed for each gene, and genes with significant differential expression are identified as marker genes for that particular bin. By applying the wilcox test to each bin along the cell-state trajectory, MeduSA can identify marker genes that are specific to each stage of the trajectory.

- gam-wald test: MeduSA associates genes along the cell-state trajectory using the generalized additive model (gam). Only genes with a false discovery rate (FDR) adjusted p-value less than 0.01 are considered. These significant genes are then ranked based on their association strength, which allows for the identification of the most relevant genes that are associated with the cell-state trajectory. To prevent certain cell-states from being overrepresented, MeduSA divides the cell-state trajectory into a specified number of intervals. Each gene is then assigned to the interval in which it has the highest mean expression. For each interval, a set of top informative genes is selected as signature genes.

To use these two methods, users can specificy the `method` in the `MeDuSA` function as `wilcox` or `gam`. Alternatively, users can use the function of `MeDuSA_marker` to select the marker genes before running cell-state deconvolution analysis. 

```r
##Set the gene selection method in MeDuSA function 
library(MeDuSA)
#wilcox
MeDuSA_obj = MeDuSA(bulk,sce,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
#gam-wald
MeDuSA_obj = MeDuSA(bulk,sce,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "gam",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)
		    
##Select marker genes using MeDuSA_marker before running deconvolution analysis
marker = MeDuSA_marker(sce[,which(sce$cell_type=='mon')],bulk,
                               GeneNumber = 200,nbins = 10,
                               family ="gaussian",k = 10,ncpu = 2,method = "wilcox")
##Documentations			       
help(MeDuSA_marker)			       
```
### 3. How to include other cell types as covariates
To account for potential confounding factors casued by other cell types, MeDuSA allows (recommends) users to include them as covariates. 

```r

```

### 4. How to use the mode of conditional auto-regressive (CAR)

### 5. Do I need to normalize the data 
Yes, we recommend the user to normalize their reference scRNA-seq data and bulk RNA-seq data into the same scale before running deconvolution analysis. 

### 6. How to get the p-value of the random-effects component


## Prepare reference data

## Compare the estimated cell-state abundance to the expected truth



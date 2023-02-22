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
bulk = readRDS("../Monocytes_bulk.rds")
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

sce@assays$RNA@counts[1:3,1:3]
              A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1
FO538757.2                             .                          .                          .
AP006222.2                             .                          .                          .
RP4-669L17.10                          .                          .                          .

sce$cell_trajectory[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                 0.9220488                  0.5167408                  0.4567616 
		 
sce$cell_type[1:3]
A/A.rds_AAACCTGCAGCGAACA-1 A/A.rds_AAACCTGGTCGACTGC-1 A/A.rds_AAACCTGGTCGCTTCT-1 
                     "mon"                      "mon"                      "mon" 
```
For compatibility with MeDuSA, the reference scRNA-seq data must be in the Seurat object format. Specifically, the reference data should be stored in `sce@assays$RNA@counts`, the cell-state trajectory in `sce$cell_trajectory`, and the cell-type in `sce$cell_type`. For more information about Seurat, please refer to the following [resource](https://satijalab.org/seurat/).



## Cell-state deconvolution analysis
```r
library(MeDuSA)

#Documentations
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
MeDuSA offers users the flexibility to input their own cultivated marker genes. In addition, MeDuSA provides two methods for selecting marker genes that are representative of the cell-state trajectory.

- wilcox test: MeDuSA divides the cells in the trajectory into a specified number of bins and applies the wilcox test, which is implemented in the `Seurat::FindMarkers` function, to each bin. The wilcox test is used to compare gene expression levels between the cells in the current bin and all other cells in the trajectory. Genes with significant differential expression are identified as marker genes for that particular bin. By performing the wilcox test for each gene at every bin along the cell-state trajectory, MeduSA can identify marker genes that are specific to each stage of the trajectory.

- gam-wald test: MeduSA uses the generalized additive model (gam) to associate genes along the cell-state trajectory, and considers only genes with an FDR-adjusted p-value less than 0.01. These significant genes are ranked based on their association strength, allowing for the identification of the most relevant genes that are associated with the cell-state trajectory. To prevent certain cell-states from being overrepresented, MeduSA divides the cell-state trajectory into a specified number of intervals and assigns each gene to the interval in which it has the highest mean expression. For each interval, a set of top informative genes is selected as signature genes.

Users can specify the `method` in the `MeDuSA` function as either `wilcox` or `gam` to utilize these two methods. Alternatively, users can select the marker genes using the `MeDuSA_marker` function before running the deconvolution analysis. 

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
                               family ="gaussian",k = 10,ncpu = 4,method = "wilcox")
##Documentations			       
help(MeDuSA_marker)			       
```
### 3. How to include other cell types as covariates
To address the possibility of confounding factors arising from other cell types, MeDuSA allows (suggests) that users include these cell types as covariates. MeDuSA provides two ways to import the other cell types in the model. 
```r
##1)--We recommend that users build the covariates matrix before running the deconvolution analysis, as this can help to save memory during the analysis. 

#1.1 load the data
sce_otherCT = readRDS("../Monocytes_OtherCell.rds")
cov_otherCT = Seurat::AverageExpression(object = sce_otherCT,group.by = 'cell_type',assays ='RNA',slot='counts')$RNA
remove(sce_otherCT)

#1.2 To input the covariates matrix into MeDuSA, users can specify the parameter of fixCov. 
MeDuSA_obj = MeDuSA(bulk,sce,fixCov = cov_otherCT,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	



##2)--Alternatively, users can also choose to merge the data of the focal cell-type and other cell-types into a single Seurat object.

#2.1 load the data
sce_otherCT = readRDS("../Monocytes_OtherCell.rds")
sce_big = merge(sce,sce_otherCT)
cell_type = c(sce$cell_type,sce_otherCT$cell_type)
sce_big$cell_type = cell_type[colnames(sce_big)]
remove(sce_otherCT); remove(sce)

#2.2 MeDuSA will automatically construct the covariates matrix based on the cell-type labels stored in sce_big$cell_type.
MeDuSA_obj = MeDuSA(bulk,sce = sce_big,
		    select.ct = 'mon',markerGene = NULL,span = 0.35,method = "wilcox",
		    resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
```

### 4. How to use the mode of conditional autoregressive (CAR)
In many LMM applications, random effects are assumed to be independent and identically distributed. However, the abundances of cells at adjacent states are likely to be correlated. MeDuSA incorporates the CAR model in the LMM to accommodate such correlations. Users can set the parameter of `CAR` as TRUE to use the CAR mode. Further, MeDuSA provides users with the ability to set the range of possible correlation strengths through the `phi` parameter.  Based on this range, MeDuSA then automatically searches for the optimal correlation strength that maximizes the likelihood function. 
```R
#The default phi is c(0.2,0.4,0.6,0.9)
MeDuSA_obj = MeDuSA(bulk,sce, CAR = TRUE,
                   select.ct = 'mon',markerGene = NULL,span = 0.35,
		   resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	

#Change the parameter space of phi
phi = c(0,0.1,0.99)
MeDuSA_obj = MeDuSA(bulk,sce, CAR = TRUE, phi = phi,
                   select.ct = 'mon',markerGene = NULL,span = 0.35,
		   resolution = 50,smooth = TRUE,fractional = TRUE,ncpu = 4)	
``` 
It is important to note that when using the CAR mode, the computational speed can become slow. This is primarily due to the computational burden of inverting the covariance matrix between cells, which can become especially significant when using large reference datasets.

### 5. How to normalize the data 
Before running the deconvolution analysis, we recommend that users normalize the reference data and the bulk data to the same scale. It is important to note that MeDuSA <big>does not</big> perform any normalization for the input reference and bulk data due to the variety in data scale, which may include raw counts, counts per million (CPM), transcripts per million (TPM), fragments per kilobase of transcript per million (FPKM), or log-transformed data.  While MeDuSA is generally robust to different scales, the heterogeneity in data scale between the bulk and reference data may negatively impact the performance. Therefore, users must carefully check and perform the appropriate normalization of their data before running MeDuSA to ensure accurate and reliable results. For example, in this tutorial, we have normalized the data into CPM scale.
```r
###To prevent exceeding the largest upper limit in R during REML iteration, we normalized the data to 1e+3 instead of 1e+6. 
bulk = sweep(bulk,2,colSums(bulk),'/')*1000
sce@assays$RNA@counts = sweep(sce@assays$RNA@counts,2,colSums(sce@assays$RNA@counts),'/')*1000

### Users can try normalizing data to other scales as well, such as the count or log-transformed scale. 
```


### 6. How to get the p-value of the random effects component


## Prepare reference data

## Compare the estimated cell-state abundance to the expected truth



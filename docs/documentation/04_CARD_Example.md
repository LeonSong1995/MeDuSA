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
First, we jointly visualize the cell type proportion matrix through scatterpie plot. Note that here because the number of spots is relatively small, so jointly visualize the cell type proportion matrix in the scatterpie plot format is duable. We do not recommend users to visualize this plot when the number of spots is > 500. Instead, we recommend users to visualize the proportion directly, i.e., using the function CARD.visualize.prop(). Details of using this function see the next example.  
```r
## set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
    "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,spatial_location = CARD_obj@spatial_location, colors = colors)
print(p1)
```
Here is an example output: 
![Example_Pie](Example_analysis_visualizePie.png)

Then, we can select some interested cell types to visualize separately. 

```r
## select the cell type that we are interested
ct.visualize = c("Acinar_cells","Cancer_clone_A","Cancer_clone_B","Ductal_terminal_ductal_like","Ductal_CRISP3_high-centroacinar_like","Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic","Fibroblasts")
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
	proportion = CARD_obj@Proportion_CARD,        
	spatial_location = CARD_obj@spatial_location, 
	ct.visualize = ct.visualize,                 ### selected cell types to visualize
	colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
	NumCols = 4)                                 ### number of columns in the figure panel
print(p2)

```
Here is an example output: 
![Example_Prop](Example_analysis_visualizeProp.png)

### 4. Visualize the cell type proportion correlation 
```r
p3 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p3)
```
Here is an example output: 
<p align="left"> 
<img src="Example_analysis_visualizeCor.png" width="700">
</p>
## Refined spatial map
A unique feature of CARD is its ability to model the spatial correlation in cell type composition across tissue locations, thus enabling spatially informed cell type deconvolution. Modeling spatial correlation allows us to not only accurately infer the cell type composition on each spatial location, but also impute cell type compositions and gene expression levels on unmeasured tissue locations, facilitating the construction of a refined spatial tissue map with a resolution much higher than that measured in the original study.
Specifically, CARD constructed a refined spatial map through the function `CARD.imputation`. The essential inputs are:

- CARD_object: CARD Object with estimated cell type compositions on the original spatial resolved transcriptomics data.
- NumGrids: Initial number of newly grided spatial locations. The final number of newly grided spatial locations will be lower than this value since the newly grided locations outside the shape of the tissue will be filtered. 
- ineibor: Numeric, number of neighbors used in the imputation on newly grided spatial locations, default is 10.  

Briefly, CARD first outlined the shape of the tissue by applying a two-dimensional concave hull algorithm on the existing locations, then perform imputation on the newly grided spatial locations. We recommend to check the exisiting spatial locations to see if there are outliers that are seriously affect the detection of the shape. 

### 1. Imputation on the newly grided spatial locations
```r
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)
## The rownames of locations are matched ...
## Make grids on new spatial locations ...
```
The results are store in `CARD_obj@refined_prop` and `CARD_obj@refined_expression`

```r
## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
	y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
library(ggplot2)
p4 <- ggplot(location_imputation, 
       aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    legend.position="bottom",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p4)
```
Here is the example of the newly grided spatial locations. The shape outline is detected well. 
<p align="center"> 
<img src="Example_analysis_grids.png" width="400">
</p>

### 2. Visualize the cell type proportion at an enhanced resolution
Now we can use the same `CARD.visualize.prop` function to visualize the cell type proportion at the enhanced resolution. But this time, the input of the function should be the imputed cell typr propotion and corresponding newly grided spatial locations.
```r                                   
p5 <- CARD.visualize.prop(
	proportion = CARD_obj@refined_prop,                         
	spatial_location = location_imputation,            
	ct.visualize = ct.visualize,                    
	colors = c("lightblue","lightyellow","red"),    
	NumCols = 4)                                  
print(p5)
```
![Example_grids](Example_analysis_grid_prop.png)

### 3. Visualize the marker gene expression at an enhanced resolution
After we obtained cell type proportion at the enhanced resolution by CARD, we can predict the spatial gene expression at the enhanced resolution. The following code is to visualize the marker gene expression at an enhanced resolution.
```r                                   
p6 <- CARD.visualize.gene(
	spatial_expression = CARD_obj@refined_expression,
	spatial_location = location_imputation,
	gene.visualize = c("Tm4sf1","S100a4","Tff3","Apol1","Crisp3","CD248"),
	colors = NULL,
	NumCols = 6)
print(p6)
```
![Example_grid_Gene](Example_analysis_grid_Gene.png)

Now, compare with the original resolution, CARD facilitates the construction of a refined spatial tissue map with a resolution much higher than that measured in the original study.
```r                                   
p7 <- CARD.visualize.gene(
	spatial_expression = CARD_obj@spatial_countMat,
	spatial_location = CARD_obj@spatial_location,
	gene.visualize = c("Tm4sf1","S100a4","Tff3","Apol1","Crisp3","CD248"),
	colors = NULL,
	NumCols = 6)
print(p7)
```
![Example_grid_Gene](Example_analysis_Original_Gene.png)

## Extension of CARD in a reference-free version: CARDfree
We extended CARD to enable reference-free cell type deconvolution and eliminate the need for the single-cell reference data. We refer to this extension as the reference-free version of CARD, or simply as CARDfree. Different from CARD, CARDfree no longer requires an external single-cell reference dataset and only needs users to input a list of gene names for previously known cell type markers. We use the same exmple dataset to illustrate the use of CARDfree. In addition to the exmple dataset, CARDfree also requires the input of marker gene list, which is in a list format with each element of the list being the cell type specific gene markers. The example marker list for runing the tutorial can be downloaded in this [page](https://yingma0107.github.io/CARD/documentation/03_data.html). 

Similar to CARD, we will first need to create a CARDfree object with the spatial transcriptomics dataset and the marker gene list 

### 1. Create an `CARDfree` object
The CARDfree object is created by the function `createCARDfreeObject`. Briefly, the essential inputs are the same as the function `createCARDObject`, except that this function does not require the single cell count and meta information matrix. Instead, it requires a markerList. 

```r
## load the marker gene list
load("./markerList.RData")
CARDfree_obj = createCARDfreeObject(
	markerList = markerList,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	minCountGene = 100,
	minCountSpot =5) 
```
The spatial data are stored in `CARDfree_obj@spatial_countMat` and `CARDfree_obj@spatial_location` while the marker list is stored in CARDfree_objj@markerList in the format of list. 
### 2. Deconvolution using CARDfree
```r
## deconvolution using CARDfree
CARDfree_obj = CARD_refFree(CARDfree_obj)
```
The results are stored in `CARDfree_obj@Proportion_CARD`. 
```r
## One limitation of reference-free version of CARD is that the cell types inferred from CARDfree do not come with a cell type label. It might be difficult to interpret the results. 
print(CARDfree_obj@Proportion_CARD[1:2,])
            CT1          CT2          CT3          CT4          CT5
10x10 0.1293363 8.155678e-12 1.342871e-10 7.934516e-06 1.029232e-01
10x13 0.7580073 1.241726e-23 1.981989e-27 3.367648e-40 1.577430e-37
               CT6          CT7          CT8          CT9         CT10
10x10 6.081635e-02 3.932038e-04 4.996482e-07 1.862520e-08 4.092203e-12
10x13 2.162971e-27 1.934013e-25 3.096247e-87 1.903368e-08 2.520620e-02
            CT11       CT12       CT13        CT14         CT15         CT16
10x10 0.05317422 0.22235226 0.05971717 0.174335930 5.853567e-38 2.537400e-10
10x13 0.02868479 0.08371092 0.09967332 0.004615396 4.214875e-37 1.017432e-04
              CT17         CT18         CT19         CT20
10x10 8.433169e-11 1.219174e-01 1.077052e-16 7.502541e-02
10x13 3.462751e-07 7.941545e-64 9.548596e-20 3.760342e-13
```

### 3. Visualization of the results of CARDfree
Note that here because the number of spots is relatively small, so jointly visualize the cell type proportion matrix in the scatterpie plot format is duable. We do not recommend users to visualize this plot when the number of spots is > 500. Instead, we recommend users to visualize the proportion directly, i.e., using the function CARD.visualize.prop(). 
```r
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
### In order to maximumply match with the original results of CARD, we order the colors to generally match with the results infered by CARD
CARDfree_obj@Proportion_CARD = CARDfree_obj@Proportion_CARD[,c(8,10,14,2,1,6,12,18,7,13,20,19,16,17,11,15,4,9,3,5)]
colnames(CARDfree_obj@Proportion_CARD) = paste0("CT",1:20)
p8 <- CARD.visualize.pie(CARDfree_obj@Proportion_CARD,CARDfree_obj@spatial_location,colors = colors)
print(p8)
```
![Example_CARDfree](Example_analysis_CARDfree_prop.png)
## Extension of CARD for single cell resolution mapping
We also extended CARD to facilitate the construction of single-cell resolution spatial transcriptomics from non-single-cell resolution spatial transcriptomics. Details of the algorithm see the main text. Briefly, we infer the single cell resolution gene expression for each measured spatial location from the non-single cell resolution spatial transcriptomics data based on reference scRNaseq data we used for deconvolution. The procedure is implemented in the function `CARD_SCMapping`. The essential inputs are:
- CARD_object: CARD object create by the createCARDObject function. This one should be the one after we finish the deconvolution procedure
- shapeSpot: a character indicating whether the sampled spatial coordinates for single cells locating in a Square-like region or a Circle-like region. The center of this region is the measured spatial location in the non-single cell resolution spatial transcriptomics data. The default is "Square", and the other option is "Circle"
- numCell: a numeric value indicating the number of cores used to accelerate the procedure.
```r
#### Note that here the shapeSpot is the user defined variable which indicates the capturing area of single cells. Details see above.
scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Square",numCell=20,ncore=10)
print(scMapping)
### the output
class: SingleCellExperiment 
dim: 16381 8560 
metadata(0):
assays(1): counts
rownames(16381): A1BG A1CF ... ZZEF1 ZZZ3
rowData names(1): rownames(count_CT)
colnames(8560): Cell686:10x10:9.97518192091957x9.83765210071579
  Cell734:10x10:10.1896061601583x9.94081321195699 ...
  Cell443:9x33:9.4747460691724x32.5644472888671
  Cell1488:9x33:9.43348842463456x33.4998327996582
colData names(7): x y ... CT Cell
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
### spatial location info and expression count of the single cell resolution data
library(SingleCellExperiment)
MapCellCords = as.data.frame(colData(scMapping))
count_SC = assays(scMapping)$counts
```
The results are stored in a SingleCellExperiment object with mapped single cell resolution counts stored in the assays slot and the information of the spatial location for each single cell as well as their relashionship to the original measured spatial location is stored in the colData slot. 

Next, we visualize the cell type for each single cell with their spatial location information 
```r
library(ggplot2)
df = MapCellCords
colors = c("#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
	"#D5CFD6","#86B1CD","#CEB28B","#EDBC63","#C59CC5","#C09CBF","#C2D567","#C9DAC3","#E1EBA0",
	"#FFED6F","#CDD796","#F8CDDE")
p9 = ggplot(df, aes(x = x, y = y, colour = CT)) + 
    geom_point(size = 3.0) +
    scale_colour_manual(values =  colors) +
    #facet_wrap(~Method,ncol = 2,nrow = 3) + 
        theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
              panel.background = element_rect(colour = "white", fill="white"),
              plot.background = element_rect(colour = "white", fill="white"),
    legend.position="bottom",
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 13,face="bold"),
    legend.text=element_text(size = 12),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm'),
    strip.text = element_text(size = 15,face="bold"))+
                                guides(color=guide_legend(title="Cell Type"))
print(p9)
```
![Example_grids](Example_analysis_scMapping.png)

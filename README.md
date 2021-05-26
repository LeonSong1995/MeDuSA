# MCTD: Mixed model-based Cell Trajectory Deconvolution. 
**Author: Liyang Song <songliyang@westlake.edu.cn>**    


## Description
MCTD is a fine-resolution deconvolution method used to predict cell abundance along a predefined cell trajectory in the given bulk RNA-seq data.



## Installation
```R
install.packages("devtools")
devtools::install_github("LeonSong1995/MCTD", build_vignettes=F)
```



## Usage

The function for **cell-type level** deconvolution in this package is `CTdcv`. It needs:  
1. Bulk RNA-Seq data. Matrix of the gene expression (count/rpkm/tpm) from the samples for which to estimate cell-type
proportions;  
2. Single-cell RNA-Seq data. [Seurat](https://satijalab.org/seurat/) object (count/rpkm/tpm) of the reference single-cell RNA-Seq data;  
3. Signature genes. We summarized signature genes for human 64 cell-types.  
**Note**: When inputting rpkm/tpm matrix (without cell size inputted), MLM can only estimate relative cell-type proportions which are not comparable among different cell-types.  

The function for **single-cell level** deconvolution in this package is `SCdcv`. It needs:  
1. Bulk RNA-Seq data. Matrix of the gene expression (count/rpkm/tpm) from the samples for which to estimate cell-type
proportions;  
2. Single-cell RNA-Seq data. [Seurat](https://satijalab.org/seurat/) object (count/rpkm/tpm) of the reference single-cell RNA-Seq data;  



## Example
MLM attached example data to show how to use: 
```R
library(MLM)
#When you load MLM, you will get these attached example data:
#1.example.bulk: human pancreas bulk gene expression data (18 samples, count-matrix);
#2.example.sce: human pancreas scRNA-Seq data (Seurat object, count-matrix);
#3.example.gene: human pancreas signatue genes (284 genes);

#basic usage:
#for cell-type level deconvolution
out.ct = CTdcv(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type = 'count')

#for single-cell level deconvolution 
out.sc = SCdcv(bulk = example.bulk,sce = example.sce,select.ct = 'beta')

```


## Contact
If you have any technical or other issues during the useage, please contact <songliyang@westlake.edu.cn>.



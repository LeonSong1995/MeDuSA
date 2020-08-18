# MLM
**Author: Liyang Song <liyang.song@ifar.ac.cn>**    
**Advisor: Jian Yang and Xiwei Sun**    
**Date: 2020-8-18**

## Description
Package implementing the ***mixed linear model (MLM)*** to estimate cell-type proportions in bulk gene expression data.
It is based on reference gene expression from single-cell RNA-Seq data. 



## Installation
```R
install.packages("devtools")
devtools::install_github("LeonSong1995/MLM", build_vignettes=TRUE)

# MLM depends on packages of "SingleCellExperiment"" and "sva"
# If you did not install these two packages,please install and library them before using MLM:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("sva")

library(SingleCellExperiment)
library(sva)
```



## Usage
The main function in this package is `Estimate`. It needs:  
1. Matrix of the gene expression (count/rpkm/tpm) from the samples for which to estimate cell-type
proportions;  
2. [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object (count/rpkm/tpm) of the reference single-cell RNA-Seq data;  
3. Signature genes.  
**Note**: When inputting rpkm/tpm matrix (without cell size inputted), MLM can only estimate relative cell-type proportions which are not comparable among different cell-types.  
```R
# library(MLM) ## If the package isn't loaded (or use MLM::Estimate and so on).
out = Estimate(bulk = bulk,sce = sce,gene = SignatureGene,data_type = 'count')
```
`out` is a list containing cell-type proportions, corresponding p-value (based on χ²(df=1))and cell-type size in each samples.


## Example
MLM attached example data to show how to use: 
```R
library(MLM)
#When you load MLM, you will get these attached example data:
#1.example.bulk: human pancreas bulk gene expression data (18 samples, count-matrix);
#2.example.sce: human pancreas scRNA-Seq data (SingleCellExperiment object, count-matrix);
#3.example.gene: human pancreas signatue genes (284 genes);
#4.ct.real: real cell-type proportions in the samples of example.bulk. 

#basic usage:
out = Estimate(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type = 'count')
out$ct.pro
out$ct.pro.p
out$cellSize
#compare the estimated cell-type proportions with real cell-type proportions
cor(ct.es$ct.pro[,colnames(ct.real)],ct.real)
#If everything goes well, you will get: 
          alpha        beta       delta       gamma
alpha  0.8920135 -0.80780295 -0.01491614 -0.41623008
beta  -0.8254393  0.90005085 -0.03776426  0.17128660
delta -0.1422445 -0.04725818  0.82439223  0.08491228
gamma -0.3996521 -0.06272664 -0.19399841  0.88499310
```


## Multiple Random Components
MLM supports fitting the mixed linear model with multiple random components.You can reset 'RanSplit' to make the model include multiple random components. 
```R
#Run model with multiple random components based on different cell-type. 
out = Estimate(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type='count',RanSplit = 'cellType')
```


## Run Parallelly
MLM supports multithreaded estimation.You can reset 'DoParallel' and 'ncpu' to let the model run parallelly.
```R
#run with 4 threads:
out = Estimate(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type='count',DoParallel = T,ncpu = 4)

#Note: If you run the model parallelly, the warning information will not be displayed on the console. You will get a RunLog file recording the running information under your current working directory. 

"~/RunLog/time_pid_logger"

```



## Warning Information
You may get some feedback warnings when the model does not work well.
```R
"1. Log-likelihood not converged. The results are not reliable. You can specify the parameter max_iter to allow for more iterations."
# If you get this warning, your model was not converged and you need to add more iterations to make the model converge.
out = Estimate(bulk = example.bulk,sce = example.sce,gene = example.gene,data_type='count',max_iter=1e+4)

"2. The V matrix is not positive. Switch to multiple linear regression!"
# If you get this warning, MLM can not work for the data you inputted and the cell-type proportions will be estimated by multiple linear regression, which may not as accurate as MLM. 

"3. X^t * V^-1 * X matrix is not invertible. Switch to multiple linear regression!"
# If you get this warning, MLM can not work for the data you inputted and the cell-type proportions will be estimated by multiple linear regression, which may not as accurate as MLM. 
```

## Contact
If you have any technical or other issue during the useage, please contact <liyang.song@ifar.ac.cn>.




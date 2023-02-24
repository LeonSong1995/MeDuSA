---
layout: page
title: Installation
description: ~
---

`MeDuSA` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 3.5.0
* R packages: Rcpp, foreach, Seurat, doSNOW, mgcv, RcppEigen, parallel, stats, BiocParallel

#### 1. Installing Dependent Packages
```r
install.packages("devtools")
install.packages("Seurat")
BiocManager::install("BiocParallel")
```

#### 2. Installing `MeDuSA`
```r
devtools::install_github("LeonSong1995/MeDuSA", build_vignettes=F)
```
#### 3. Loading package
```r
library(MeDuSA)
```

This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- Windows 10: Home (1903)
- MAC: OSX (11.2)
- Linux: CentOS-8

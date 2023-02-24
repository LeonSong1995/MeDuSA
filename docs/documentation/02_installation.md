---
layout: page
title: Installation
description: ~
---

`CARD` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.0.0.
* R packages: SingleCellExperiment, SummarizedExperiment, concaveman, sp, Matrix, methods, ggplot2, ggcorrplot, MuSiC, fields, MCMCpack, dplyr, sf, RANN, stats, reshape2, RColorBrewe, scatterpie, grDevices, stats, nnls, pbmcapply, spatstat, gtools, RcppML, NMF


#### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

#### 2. Install `MeDuSA`
```r
devtools::install_github("LeonSong1995/MeDuSA", build_vignettes=F)
```
#### 3. Load package
```r
library(MeDuSA)
```

This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- Windows 10: Home (1903)
- MAC: OSX (10.14.1)
- Linux: Ubuntu (16.04.6)

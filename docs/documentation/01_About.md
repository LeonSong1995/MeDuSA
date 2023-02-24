---
layout: page
title: About
---

Parameters of `MeDuSA`
-------------------
\code{bulk}: A matrix of bulk RNA-seq data. Each row corresponds to a specific gene, and each column corresponds to a particular sample.
\code{sce}: A \code{seurat} object of the reference scRNA-seq data (see \code{\link{seurat}}). Meta-data of the seurat object (\code{sce@meta.data}) needs to include two columns named as \code{cell_type} and \code{cell_trajectory}.
\code{select.ct}: A character variable of the target cell type.
\code{ncpu}: The number of CPU cores to be used.
\code{smooth}: A Boolean variable to determine whether to smooth the predicted cell-state abundance or not. The default value is TRUE.
\code{smoothMethod}: A character variable to determine the smoothing method. The default method is \code{\link{loess}}.
\code{nbins}: A numeric variable to determine the number of bins along the cell trajectory, which is used to ensure the selected genes are uniformly scattered along the given trajectory. The default value is 10.
\code{resolution}: A numeric variable to determine the number of cell bins along the cell trajectory (i.e., the resolution of the deconvolution analysis). The default value is 50.
\code{knots}: A numeric variable to specify the number of knots. The default value is 10.
\code{start}: A numeric vector for the initial value of the restricted maximum likelihood (REML). The default value is c(1e-5,1e-2).
\code{family}: A character variable to specify the distribution of GAM. See \code{\link{family.mgcv}} for a full list of what is available.
\code{maxiter}: The iteration number of REML. The default value is 1e+4.
\code{adj}: A Boolean variable to determine whether to include covariates when predicting the cell-state abundance.
\code{CAR}: A Boolean variable to determine whether to model abundance correlations among cells.
\code{phi}: A numeric vector for searching the optimal cell correlations. The default value is c(0.2,0.4,0.6,0.9).
\code{fixCov}: A matrix (vector) of fixed covariates in the model (i.e., covariates for estimating cell-state abundance). The default value is NULL.
\code{markerGene}: A character vector of marker genes over the cell-state trajectory. The default value is NULL. With the default, \code{MeDuSA} selects genes using the GAM-Wald or Wilcoxon rank-sum test.
\code{method}: A character variable to specify the method used in selecting marker genes (Wilcoxon or GAM). The default value is Wilcoxon.
\code{GeneNumber}: A numeric variable to determine the number of signature genes. The default value is 200.
\code{span}: A numeric variable to control the degree of smoothing for \code{\link{loess}}.
\code{neighbor}: A numeric variable to determine the number of neighboring cells used in smoothing (do not use when the smooth method is loess).
\code{fractional}: A Boolean variable to determine whether to normalize the estimated cell-state abundance to the fractional abundance (0,1).


Cite `MeDuSA`
-------------------
Liyang Song, Xiwei Sun, Ting Qi, Jian Yang

Contact
-------------------
If you have any questions about MeDuSA, please feel free to leave messages on the github [issues](https://github.com/LeonSong1995/MeDuSA/issues) or contact me: songliyang@westlake.edu.cn.

Our group
-------------------
[Yang lab website](https://yanglab.westlake.edu.cn/)


<!-- README.md is generated from README.Rmd. Please edit that file -->

# CLUEY

<!-- badges: start -->
<!-- badges: end -->

This is an R package for estimating the number of clusters in uni and
multi-modal single-cell data. CLUEY uses cell-type identity markers to
guide the clustering process and performs recursive clusters to ensure
that sub-populations are captured.

| ![](man/figures/CLUEY_figure.png) |
|:---------------------------------:|

## Dependencies

CLUEY requires both keras and tensorflow, please have both installed.
You can follow the instructions provided at this
[link](https://tensorflow.rstudio.com/install/).

## Installation

CLUEY can be installed using the following command:

``` r
library(devtools)
install_github("SydneyBioX/CLUEY")
```

## Generating knowledge base

You can generate your own knowledge base using the
`generate_knowledgebase` function like below:

``` r
knowledgeBase <- generate_knowledgebase(reference=logcounts(sce), celltypes=sce$cellType)
```

## Cluster data

In this example, we will upload an example knowledge base generated from
the [Mouse Cell Atlas
(FACS)](https://www.nature.com/articles/s41586-020-2496-1) and cluster
an example query dataset which was subsampled from [Zilionis et
al.](https://pubmed.ncbi.nlm.nih.gov/30979687/) using the `run_CLUEY`
function.

``` r
library(CLUEY)
library(scater)
library(ggplot2)
library(gridExtra)
set.seed(3435)

# Load example knowledge base
data(mcaFACS)

# Load example query data
data(exampleData)

# Run CLUEY
# If your logcounts matrix is in dgCMatrix format, then you'll need to convert it to a matrix using `as.matrix()`
clustering_results <- run_CLUEY(rna=as.matrix(logcounts(exampleData)), knowledgebase=mcaFACS, k_limit=10)
#> 50/50 - 0s - 73ms/epoch - 1ms/step
#> 6/6 - 0s - 41ms/epoch - 7ms/step
#> 7/7 - 0s - 45ms/epoch - 6ms/step
#> 7/7 - 0s - 50ms/epoch - 7ms/step
#> 6/6 - 0s - 43ms/epoch - 7ms/step
```

## Viewing results

We can now view the results of the clustering performed by CLUEY. CLUEY
predicts there to be 5 clusters in the data.

``` r
set.seed(3435)

# View the optimal number of clusters predicted by CLUEY
clustering_results$optimal_K
#> [1] 5

# We can store the results in the metadata of our SingleCellExperiment object. 
colData(exampleData) <- cbind(colData(exampleData), clustering_results$predictions)

# Run UMAP to visualise clusters
exampleData <- runPCA(exampleData)
exampleData <- runUMAP(exampleData)
umap <- data.frame(reducedDim(exampleData, "UMAP"))
umap$cluster <- as.factor(exampleData$cluster)
umap$correlation <- exampleData$correlation

ggplot(umap, aes(x=UMAP1, y=UMAP2, color=cluster)) + geom_point() + theme_classic()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

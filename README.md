
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CLUEY

<!-- badges: start -->
<!-- badges: end -->

This is an R package for estimating the number of clusters in uni and
multi-modal single-cell data. CLUEY uses cell-type identity markers to
guide the clustering process and performs recursive clusters to ensure
that sub-populations are captured.

## Installation

CLUEY can be installed using the following command:

``` r
library(devtools)
install_github("SydneyBioX/CLUEY")
```

## Generating reference

We provide four references but you can generate your own references
using the `generateReference` function like below:

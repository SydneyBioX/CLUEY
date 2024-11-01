% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCLUEY.R
\name{runCLUEY}
\alias{runCLUEY}
\title{runCLUEY}
\usage{
runCLUEY(
  exprsMatRNA,
  exprsMatOther = NULL,
  knowledgeBase,
  corMethod = "spearman",
  propGenes = 0.25,
  kLimit = 20,
  subK = 3,
  minCells = 20,
  recursive = TRUE,
  encodingDim1 = 50,
  encodingDim2 = 10,
  hiddenDimsMultimodal = 50,
  nEpochs = 30
)
}
\arguments{
\item{exprsMatRNA}{A numeric matrix of normalised gene expression data where rows are
features and columns are cells.}

\item{exprsMatOther}{A list of additional modalities, such as ATAC-seq or CITE-seq,
where rows are the features and columns are cells.}

\item{knowledgeBase}{knowledgeBase generated from the function `generateKnowledgeBase`.}

\item{corMethod}{Correlation method used when determining the most similar cell-types
for each cluster. Options are "Spearman" or "Pearson". Default is Spearman correlation.}

\item{propGenes}{Proportion of genes to use when calculating the correlation of the expression
profiles between each cluster and knowledgeBase cell-type. Default is 0.25.}

\item{kLimit}{The maximum number of clusters to estimate. Default is 20.}

\item{subK}{The maximum number of clusters to estimate when performing iterative clustering. Default is 3.}

\item{minCells}{The minimum number of cells that a cluster should contain. Default is 20 cells.}

\item{recursive}{A boolean specifying whether to perform clustering recursively. Default is TRUE.}

\item{encodingDim1}{Dimension of the hidden layer for autoencoder at the first level of cluster estimation. Deault is 50.}

\item{encodingDim2}{Dimension of the hidden layer for autoencoder at the second level of cluster estimation. Deault is 20.}

\item{hiddenDimsMultimodal}{Dimension of the hidden layer for multi-modal autoencoder. Ignored if `exprsMatOther` is NULL. Default is 50.}

\item{nEpochs}{Number of epochs when training the autoencoder.}
}
\value{
A list containing the optimal number of clusters `optimal_k` and
final clustering score `score`. Cluster assignments and correlation scores
for each cluster can be in `predictions`.
}
\description{
Perform cluster estimation of single or multi-modal single-cell data
}

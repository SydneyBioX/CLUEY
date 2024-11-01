% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateKnowledgeBase.R
\name{generateKnowledgeBase}
\alias{generateKnowledgeBase}
\title{generateKnowledgeBase}
\usage{
generateKnowledgeBase(
  exprsMat,
  celltypes,
  batch = NULL,
  method = "ds",
  minCells = 20
)
}
\arguments{
\item{exprsMat}{A numeric matrix of normalised gene expression data where rows are
features and columns are cells.}

\item{celltypes}{vector of cell-types. Must be the same length as the number of columns, i.e. cells,
passed to `exprsMat`.}

\item{batch}{vector of batch labels. Must be the same length as the number of columns, i.e. cells,
passed to `exprsMat`. Default is NULL.}

\item{method}{method used for identifying variable genes and ranking them. Default is `ds` for differentially stable genes.
To use differentially expressed genes, specify `de`.}

\item{minCells}{the minimum number of cells to use when pseudobulking. Cell-types with less than this number are ignored.
Default is 20 cells.}
}
\value{
A named list containing pseudobulked expression profiles for each cell-type provided in exprsMat.
}
\description{
Generate pseudobulked exprsMat for cluster estimation performed in runCLUEY
}

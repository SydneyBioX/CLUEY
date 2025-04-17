#' generateKnowledgeBase
#'
#' @description
#' Generate pseudobulked exprsMat for cluster estimation performed in runCLUEY
#'
#' @param exprsMat A numeric matrix of normalised gene expression data where rows are
#' features and columns are cells.
#' @param celltypes vector of cell-types. Must be the same length as the number of columns, i.e. cells,
#' passed to `exprsMat`.
#' @param batch vector of batch labels. Must be the same length as the number of columns, i.e. cells,
#' passed to `exprsMat`. Default is NULL.
#' @param method method used for identifying variable genes and ranking them. Default is `ds` for differentially stable genes.
#' To use differentially expressed genes, specify `de`.
#' @param minCells the minimum number of cells to use when pseudobulking. Cell-types with less than this number are ignored.
#' Default is 20 cells.
#' @return A named list containing pseudobulked expression profiles for each cell-type provided in exprsMat.
#' @export
#'
generate_knowledgebase <- function(reference, celltypes, batch = NULL, method = "ds", min_cells = 20){
  if(!(2*ncol(reference) == sum(ncol(reference), length(celltypes)))){
    stop("number of cells and length of celltypes vector do not match. Batch vector must also be the same length if it is being used...")
  }

  rownames(reference) <- toupper(rownames(reference))
  celltypes <- make.names(celltypes)
  idx <- rowSums(reference)
  idx <- which(idx == 0)

  if(length(idx) > 0){

    reference <- suppressWarnings(as.matrix(reference[-idx,]))

  }else{

    reference <- suppressWarnings(as.matrix(reference))

  }

  if(method == "ds"){

    if(!is.null(batch)){

      message("Identifying differentially stable genes per batch...\n")
      batches <- unique(batch)

      batch_genes <- lapply(batches, function(x){

        idx <- which(batch == x)

        n_celltypes <- unique(celltypes[idx])

        if(length(n_celltypes) > 1){

          keep_celltypes <- names(table(celltypes[idx]) >= min_cells)

          if(length(keep_celltypes) > 0){

            idx <- which(celltypes %in% keep_celltypes)
            reference_cepo <- suppressWarnings(Cepo::Cepo(reference[,idx], celltypes[idx], min_cells = min_cells, exprs_pct = 0.05))
            reference_top_genes <- Cepo::topGenes(reference_cepo, n = nrow(reference_cepo$stats), returnValues = T)
            ref_genes <- lapply(reference_top_genes, function(x){x[x > 0]})
            ref_genes

          }

        }

      })

      batch_genes <- unlist(batch_genes, recursive = FALSE)

      grouped_celltypes <- split(batch_genes, names(batch_genes))
      ref_genes <- lapply(grouped_celltypes, function(x){
        genes <- names(unlist(unname(x)))
        common_genes <- genes[which(duplicated(genes) == TRUE)]

        if(length(common_genes) > 0){

          average_values <- tapply(unlist(unname(x)), names(unlist(unname(x))), mean)
          average_values[order(average_values, decreasing = TRUE)]

        }else{

          unlist(unname(x))

        }

      })

    }else{

      message("Identifying differentially stable genes...\n")

      if(sum(table(celltypes) < min_cells) > 0){

        stop(paste0("Less than ", min_cells, " cells for at least one cell-type, please filter out or decrease mincells parameter..."))

      }

      reference_cepo <- suppressWarnings(Cepo::Cepo(reference, celltypes, min_cells = min_cells, exprs_pct = 0.05))
      reference_top_genes <- Cepo::topGenes(reference_cepo, n = nrow(reference_cepo$stats), returnValues = T)
      ref_genes <- lapply(reference_top_genes, function(x){x[x > 0]})

    }

  }

  if(method == "de"){

    message("Identifying differentially expressed genes...\n")
    seurat_ref <- suppressWarnings(CreateSeuratObject(counts = reference, assay = "RNA", meta.data = data.frame(celltypes = celltypes)))
    seurat_ref <- Seurat::NormalizeData(seurat_ref)
    Idents(object = seurat_ref) <- celltypes
    tmp <- Seurat::FindAllMarkers(seurat_ref, test.use="t")

    ref_genes <- lapply(unique(tmp$cluster), function(celltype){

      ct_tmp <- tmp[tmp$cluster %in% celltype,]
      ct_tmp <- ct_tmp[order(ct_tmp$p_val_adj, decreasing = FALSE),]
      ct_tmp[ct_tmp$p_val_adj < 0.05,]$gene

    })

    names(ref_genes) <- make.names(levels(seurat_ref))

  }

  reference_pseudobulk <- mapply(function(genes, celltype){

    if(method == "ds"){

      genes <- names(genes)

    }else{

      genes <- intersect(genes, rownames(reference))

    }

    tmp <- reference[genes, which(celltypes == celltype)]
    tmp <- rowMeans(tmp)
    tmp

  },ref_genes, names(ref_genes))



  return(reference_pseudobulk)

}

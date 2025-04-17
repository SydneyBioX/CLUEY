#' runCLUEY
#'
#' @description
#' Perform cluster estimation of single or multi-modal single-cell data
#'
#' @param rna A numeric matrix of normalised gene expression data where rows are
#' features and columns are cells.
#' @param modalities A list of additional modalities, such as ATAC-seq or CITE-seq,
#' where rows are the features and columns are cells.
#' @param knowledgebase knowledgebase generated from the function `generateKnowledgeBase`.
#' @param cor_method Correlation method used when determining the most similar cell-types
#' for each cluster. Options are "Spearman" or "Pearson". Default is Spearman correlation.
#' @param prop_genes Proportion of genes to use when calculating the correlation of the expression
#' profiles between each cluster and knowledgeBase cell-type. Default is 0.25.
#' @param k_limit The maximum number of clusters to estimate. Default is 20.
#' @param sub_k The maximum number of clusters to estimate when performing iterative clustering. Default is 3.
#' @param min_cells The minimum number of cells that a cluster should contain. Default is 20 cells.
#' @param recursive A boolean specifying whether to perform clustering recursively. Default is TRUE.
#' @param nEpochs Number of epochs when training the autoencoder.
#'
#' @return A list containing the optimal number of clusters `optimal_k` and
#' cluster assignments and correlation scores for each cluster can be in `predictions`.
#' @export
#'

run_CLUEY <-  function(rna, modalities=NULL, knowledgebase, clus_method="spectral", cor_method = "spearman", prop_genes=0.25, k_limit=20, sub_k=3, min_cells=30,
                       recursive=TRUE, n_epochs=50){
  unimodal <- TRUE
  stopifnot(min_cells >= 2*sub_k)
  rownames(rna) <- toupper(rownames(rna))

  if (!is.null(modalities)) {

    unimodal <-  FALSE
    rownames(modalities) <- toupper(rownames(modalities))

  }

  message("Selecting HVGs...\n")
  model <- scran::modelGeneVar(rna)
  hvg_genes <- scran::getTopHVGs(model, n=0.05*nrow(rna))
  hvg_data <- t(rna[hvg_genes, ])

  message("Initial estimate of k...\n")
  sil_width <- numeric(k_limit)
  scaled_exprs <- scale(hvg_data)
  umap <- uwot::umap(scaled_exprs)
  distMatrix <- distances::distances(umap)

  sil_width <- c()

  for(i in 2:k_limit){

    kmeans_result <- stats::kmeans(umap, centers=i, nstart=5)
    ss <- cluster::silhouette(kmeans_result$cluster, distMatrix)
    sil_width <- c(sil_width, mean(ss[,3]))

  }

  k_estimate <- which(sil_width==max(sil_width))+1

  if (unimodal) {

    message("Reducing dimensions...\n")
    reduced_dims <- get_encoding(list(hvg_data), hidden_dims = k_estimate*50, encoded_dim = k_estimate*10, epochs=n_epochs, batch_size=16)

  } else if (!unimodal) {

    message("Reducing dimensions...\n")
    reduced_dims <- get_encoding(list(hvg_data, t(modalities)), hidden_dims = k_estimate*50, encoded_dim = k_estimate*10, epochs=n_epochs, batch_size=16)

  }

  dist_matrix <- as.matrix(distances::distances(reduced_dims))
  affinity_matrix <- SNFtools::affinityMatrix(dist_matrix, K = 30, sigma = 0.4)

  k_annotations <- list()
  k_clusters <- list()

  message("Estimating k...\n")

  for(k in 2:k_limit){

    if(clus_method=="spectral"){

      cluster_res <- spectralClustering(as.matrix(affinity_matrix),  k)
      clusters <- cluster_res$label

    }else{

      cluster_res <- stats::kmeans(umap, centers=k, nstart=5)
      clusters <- cluster_res$cluster

    }

    k_clusters[[length(k_clusters) + 1]] <- clusters
    k_annotations[[length(k_annotations) + 1]] <- label_clusters(exprs_mtx=rna, clusters=clusters, cor_method=cor_method, knowledgebase=knowledgebase, prop_genes=prop_genes)

  }

  names(k_annotations) <- 2:k_limit
  names(k_clusters) <- 2:k_limit
  optimal_annotations <- find_optimal_k(k_annotations=k_annotations, k_clusters=k_clusters, dist_matrix=dist_matrix)

  cluey_df <- data.frame(cell_id = colnames(rna),
                         spectral_cluster = optimal_annotations$clusters,
                         annotation = gsub("[.]", " ", optimal_annotations[[1]]$annotation[optimal_annotations$clusters]),
                         correlation = optimal_annotations[[1]]$correlation[optimal_annotations$clusters])

  cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
  cluey_df <- cluey_df[match(colnames(rna), cluey_df$cell_id),]
  max_cluster <- max(cluey_df$cluster)

  current_result <- list(optimal_K = max_cluster, annotations = cluey_df)
  tmp_current_result <- current_result
  results <- list()

  message("Running multiscale clustering...\n")

  for(i in 1:max_cluster){

    idx <- which(cluey_df$cluster == i)
    cells <- cluey_df$cell_id[idx]

    if(length(cells) > min_cells){

      tmp_current_result[["annotations"]] <- tmp_current_result$annotations[tmp_current_result$annotations$cell_id %in% cells,]

      if (!unimodal) {

        tmp_modalities <- modalities[, cells]

      }else{

        tmp_modalities <- NULL

      }

      tmp_results <- rerun_CLUEY(rna= rna[,cells], modalities=tmp_modalities, reduced_dims = reduced_dims[cells,],
                                 knowledgebase = knowledgebase, clus_method = clus_method, cor_method = cor_method, k_limit = sub_k,
                                 previous_result = tmp_current_result, prop_genes = prop_genes, recursive = recursive,
                                 min_cells = min_cells, unimodal=unimodal, n_epochs=n_epochs)

      results[[length(results)+1]] <- tmp_results

    }

  }


  if(length(results) > 0){

    cluey_df <- current_result$annotations

    for(i in 1:length(results)){

      result <- results[[i]]

      if(!is.null(result)){

        cells <- result$annotations$cell_id
        cluey_df <- cluey_df[!(cluey_df$cell_id %in% cells),]
        cluey_df <- rbind(cluey_df, results[[i]]$annotations)

      }

    }

    cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
    cluey_df <- cluey_df[match(colnames(rna), cluey_df$cell_id),]

    # score <- FisherZInv(mean(FisherZ(unique(cluey_df$correlation))))
    # pvalue <- pnorm(mean(FisherZ(unique(cluey_df$correlation))), lower.tail = F)
    return(list(optimal_K = length(unique(cluey_df$annotation)), annotations = cluey_df))

  }

  message("Returning results...")

  return(current_result)

}

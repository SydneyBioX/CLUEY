###################################################################
# Function called during recursive clustering stage
###################################################################

rerun_CLUEY <-  function(rna, modalities=modalities, reduced_dims, knowledgebase, clus_method, cor_method, k_limit, previous_result, prop_genes, recursive, min_cells,
                         unimodal=unimodal, n_epochs=n_epochs){

  encoded_dim <- 5

  if(ncol(rna) > min_cells){

    if(ncol(reduced_dims) > encoded_dim){

      model <- scran::modelGeneVar(rna)
      hvg_genes <- scran::getTopHVGs(model, n=0.05*nrow(rna))
      hvg_data <- t(rna[hvg_genes, ])

      if (unimodal) {

        reduced_dims_tmp <- cbind(hvg_data, reduced_dims)
        reduced_dims <- get_encoding(list(reduced_dims_tmp), hidden_dims = nrow(reduced_dims_tmp)/2, encoded_dim = encoded_dim, epochs=n_epochs, batch_size = 32)

      } else if (!unimodal) {

        reduced_dims <- get_encoding(list(hvg_data, t(modalities)), hidden_dims = nrow(hvg_data)/2, encoded_dim = encoded_dim, epochs=n_epochs, batch_size =32)

      }

    }

    dist_matrix <- as.matrix(dist(reduced_dims))
    affinity_matrix <- SNFtool::affinityMatrix(dist_matrix, K = 20, sigma = 0.4)

    k_annotations <- list()
    k_clusters <- list()

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
    max_cluster <- max(cluey_df$cluster)

    current_result <- list(optimal_K = max_cluster, annotations=cluey_df)
    tmp_current_result <- current_result
    results <- list()

    if(recursive == TRUE){

      if(max_cluster > 1){

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

            tmp_results <- rerun_CLUEY(rna = rna[, cells], modalities = tmp_modalities, reduced_dims = reduced_dims[cells,],
                                       knowledgebase = knowledgebase, clus_method = clus_method, cor_method = cor_method, k_limit = k_limit,
                                       previous_result = tmp_current_result, prop_genes = prop_genes, recursive = recursive,
                                       min_cells = min_cells, unimodal=unimodal, n_epochs=n_epochs)

            results[[length(results)+1]] <- tmp_results

          }

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

        return(list(optimal_K = length(unique(cluey_df$annotation)), annotations = cluey_df))

      }

    }else{

      message("...")
      return(current_result)

    }

  }

  message("...")
  return(previous_result)

}

############################################################################
# Function for labeling cluster to most similar cell-type in knowledge base
############################################################################

label_clusters <- function(exprs_mtx, clusters, cor_method, reference, prop_genes){

  # generate pseudo-bulk for query according to clusters
  pseudobulk_profiles <- generate_pseudobulk(exprs_mtx, clusters)

  k_list <- vector("list", length = length(pseudobulk_profiles))

  for(x in 1:length(pseudobulk_profiles)){

    tmp_scores <- vector("numeric", length = length(reference))
    current_cluster <- pseudobulk_profiles[[x]]
    cluster_ID <- names(pseudobulk_profiles)[[x]]

    for(y in 1:length(reference)){

      nGenes <- length(reference[[y]])*prop_genes
      common_genes <- names(reference[[y]])[names(reference[[y]]) %in% names(current_cluster)][1:nGenes] # subset for top genes
      weight <- length(common_genes)/nGenes

      corr <- .Machine$double.xmin

      if(length(common_genes > 1)){

        corr <- cor(current_cluster[common_genes], reference[[y]][common_genes], method = cor_method)

        if(is.na(corr)){

          corr <- .Machine$double.xmin

        }else if(corr == 1){

          corr == 0.999 # can't have 1 for fisher's

        }

      }

      tmp_scores[y] <- weight*corr

    }

    names(tmp_scores) <- names(reference)
    k_list[[x]] <- tmp_scores

  }

  names(k_list) <- names(pseudobulk_profiles)

  annotations <- mapply(function(x,y){

    idx <- which.max(x) # get maximum correlation
    return(list(annotation = names(x)[idx], score = x[[idx]]))

  }, k_list, names(k_list), SIMPLIFY = FALSE)

  return(annotations)

}

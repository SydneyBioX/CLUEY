####################################################################
# Function for labeling cluster to most similar reference cell-type
####################################################################

annotateClusters <- function(exprsMat, clusters, corMethod, reference, propGenes){

  # generate pseudo-bulk for query according to clusters
  pseudobulk_profiles <- generate_pseudobulk(exprsMat, clusters)

  k_list <- vector("list", length = length(pseudobulk_profiles))

  for(x in 1:length(pseudobulk_profiles)){

    tmp_scores <- vector("numeric", length = length(reference))
    current_cluster <- pseudobulk_profiles[[x]]
    #cluster_ID <- as.numeric(gsub("X", "", names(top_exprsMat_cepo_genes)[[x]]))
    cluster_ID <- names(pseudobulk_profiles)[[x]]

    for(y in 1:length(reference)){
      nGenes <- length(reference[[y]])*propGenes
      common_genes <- names(reference[[y]])[names(reference[[y]]) %in% names(current_cluster)][1:nGenes] # subset for top genes
      weight <- length(common_genes)/nGenes
      # Calculate and pass in p-values

      corr <- .Machine$double.xmin
      if(length(common_genes > 1)){
        corr <- cor(current_cluster[common_genes], reference[[y]][common_genes], method = corMethod)
        if(is.na(corr)){
          corr <- .Machine$double.xmin
        }else if(corr == 1){
          corr == 0.99 # can't have 1 for fisher's
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

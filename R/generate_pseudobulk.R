###################################################################
# Function for generating pseudobulk profiles
###################################################################

generate_pseudobulk <- function(exprs_mtx, clusters){

  cluster_labels <- unique(clusters)
  num_clusters <- length(cluster_labels)
  pseudobulk_profiles <- vector("list", num_clusters)

  for(i in 1:num_clusters){

    profile <- exprs_mtx[, which(clusters == cluster_labels[i])]

    if(!is.null(ncol(profile))){

      pseudobulk_profiles[[i]] <- rowMeans(exprs_mtx[, which(clusters == cluster_labels[i])])

    }else{

      pseudobulk_profiles[[i]] <- profile

    }


  }

  names(pseudobulk_profiles) <- cluster_labels
  return(pseudobulk_profiles)

}

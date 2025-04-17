###################################################################
# Function for determining optimal K
###################################################################

find_optimal_k <- function(k_annotations, k_clusters, dist_matrix){

  k_scores <- lapply(k_annotations, function(x){ lapply(x, function(y){y$score}) })
  k_celltypes <- lapply(k_annotations, function(x){ lapply(x, function(y){y$annotation}) })
  avg_corrs <- vector("numeric", length = length(k_scores))
  clustering_scores <- vector("numeric", length = length(k_scores))
  ar2 <- vector("numeric", length = length(k_scores))
  n <- length(k_clusters[[1]])

  for(i in 1:length(k_scores)){

    tmp <- k_annotations[[i]]
    clusters <- k_clusters[[i]]
    k <- length(unique(clusters))
    ct_labels_tmp <- unlist(lapply(tmp, function(x){x$annotation}))
    correlations <- unlist(lapply(tmp, function(x){x$score}))

    ct_labels <- gsub("[.]", " ", ct_labels_tmp[clusters])

    avg_corr <- DescTools::FisherZ(unlist(correlations))
    avg_corr <- DescTools::FisherZInv(mean(avg_corr))
    avg_corrs[i] <- avg_corr

    if(length(unique(ct_labels))>0){

      res <- fpc::cluster.stats(d = dist_matrix, clustering = clusters)
      r_squared <- res$average.between/(res$average.within+res$average.between)

    }else{

      r_squared <- 1

    }

    adj_rsquared <- 1-(((1-r_squared)*(n-1))/(n-k^2))
    ar2[i] <- adj_rsquared
    clustering_scores[i] <- adj_rsquared*avg_corr
  }

  idx <- which.max(clustering_scores)
  optimal_annotations_tmp <- k_annotations[[idx]]
  clusters <- names(optimal_annotations_tmp)

  optimal_annotations <- unlist(lapply(optimal_annotations_tmp, function(x){x$annotation}))
  optimal_correlations <- unlist(lapply(optimal_annotations_tmp, function(x){x$score}))
  names(optimal_annotations) <- clusters
  names(optimal_correlations) <- clusters
  df <- data.frame(annotation = optimal_annotations,
                   correlation = optimal_correlations)
  df$cluster <- as.numeric(rownames(df))
  df <- df[order(df$cluster), ]

  return(list(annotations = df, clusters = k_clusters[[idx]], score = max(clustering_scores)))

}

###################################################################
# Function for determining optimal K
###################################################################

findOptimalK <- function(k_annotations, k_clusters){
  k_scores <- lapply(k_annotations, function(x){ lapply(x, function(y){y$score}) })
  k_celltypes <- lapply(k_annotations, function(x){ lapply(x, function(y){y$annotation}) })
  clustering_scores <- vector("numeric", length = length(k_scores))
  for(x in 1:length(k_scores)){
    l_scores <- unlist(k_scores[[x]])
    l_celltypes <- unlist(k_celltypes[[x]])
    clustering_score <- mean(DescTools::FisherZ(unlist(l_scores)))
    clustering_scores[x] <- DescTools::FisherZInv(clustering_score)
  }

  idx <- which.max(clustering_scores)
  optimal_annotations_tmp <- k_annotations[[idx]]
  clusters <- names(optimal_annotations_tmp)

  optimal_annotations <- unlist(lapply(optimal_annotations_tmp, function(x){ x$annotation}))
  optimal_correlations <- unlist(lapply(optimal_annotations_tmp, function(x){ x$score}))
  names(optimal_annotations) <- clusters
  names(optimal_correlations) <- clusters
  df <- data.frame(annotation = optimal_annotations,
                   correlation = optimal_correlations)
  df$cluster <- as.numeric(rownames(df))
  df <- df[order(df$cluster), ]

  return(list(annotations = df, clusters = k_clusters[[idx]], score = max(clustering_scores)))
}

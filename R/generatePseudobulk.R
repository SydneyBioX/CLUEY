
# CLUEY with weighted average correlation coefficients - Fisher's Z method
annotateClusters <- function(exprsMat, clusters, corMethod, reference, propGenes){

  # generate pseudo-bulk for query according to clusters
  pseudobulk_profiles <- generatePseudobulk(exprsMat, clusters)

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


runCLUEY <-  function(exprsMatRNA, exprsMatOther=NULL, reference, corMethod = "spearman", propGenes=0.25, kLimit=20, subK=3, minCells=20,
                      recursive=TRUE, encodingDim1=50, encodingDim2=10, hiddenDimsMultimodal=50, nEpochs=50){
  unimodal <- TRUE
  stopifnot(minCells >= 2*subK)
  rownames(exprsMatRNA) <- toupper(rownames(exprsMatRNA))

  if (!is.null(exprsMatOther)) {
    unimodal <-  FALSE
    rownames(exprsMatOther) <- toupper(rownames(exprsMatOther))
  }

  if (unimodal) {

    message("Selecting HVGs...\n")
    model <- modelGeneVar(exprsMatRNA)
    hvg_genes <- getTopHVGs(model, n=0.05*nrow(exprsMatRNA))
    hvg_data <- exprsMatRNA[hvg_genes, ]
    message("Reducing dimensions...\n")
    reducedDims <- getEncodingMultiModal(list(hvg_data), hiddenDims = encodingDim1, encodedDim = encodingDim1, epochs=nEpochs, batchSize=32)
  } else if (!unimodal) {
    message("Selecting HVGs...\n")
    model <- modelGeneVar(exprsMatRNA)
    hvg_genes <- getTopHVGs(model, n=0.05*nrow(exprsMatRNA))
    hvg_data <- exprsMatRNA[hvg_genes, ]
    message("Reducing dimensions...\n")
    reducedDims <- getEncodingMultiModal(list(hvg_data, exprsMatOther), hiddenDims = hiddenDimsMultimodal, encodedDim = encodingDim1, epochs=nEpochs, batchSize=32)
  }


  dist_matrix <- as.matrix(distances::distances(reducedDims))
  affinity_matrix <- affinityMatrix(dist_matrix, K = 30, sigma = 0.4)

  k_annotations <- list()
  k_clusters <- list()
  message("Estimating initial k...\n")
  for(k in 2:kLimit){
    cluster_res <- spectralClustering(as.matrix(affinity_matrix),  k)
    clusters <- cluster_res$label
    k_clusters[[length(k_clusters) + 1]] <- clusters
    k_annotations[[length(k_annotations) + 1]] <- annotateClusters(exprsMatRNA, clusters, corMethod = corMethod, reference = reference, propGenes = propGenes)

  }
  names(k_annotations) <- 2:kLimit
  names(k_clusters) <- 2:kLimit
  optimal_annotations <- findOptimalK(k_annotations, k_clusters) # Pass in dis matrix

  cluey_df <- data.frame(cell_id = colnames(exprsMatRNA),
                         spectral_cluster = optimal_annotations$clusters,
                         annotation = gsub("[.]", " ", optimal_annotations[[1]]$annotation[optimal_annotations$clusters]),
                         correlation = optimal_annotations[[1]]$correlation[optimal_annotations$clusters])

  cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
  cluey_df <- cluey_df[match(colnames(exprsMatRNA), cluey_df$cell_id),]
  max_cluster <- max(cluey_df$cluster)
  score <- FisherZInv(mean(FisherZ(unique(cluey_df$correlation))))
  pvalue <- pnorm(mean(FisherZ(unique(cluey_df$correlation))), lower.tail = F)

  current_result <- list(optimal_K = max_cluster, score = round(optimal_annotations$score, digits = 2), annotations = cluey_df, pvalue = round(pvalue, digits = 2))
  tmp_current_result <- current_result
  results <- list()

  message("Running multiscale clustering...\n")
  for(i in 1:max_cluster){

    idx <- which(cluey_df$cluster == i)
    cells <- cluey_df$cell_id[idx]
    if(length(cells) > minCells){

      tmp_current_result[["annotations"]] <- tmp_current_result$annotations[tmp_current_result$annotations$cell_id %in% cells,]

      if (!unimodal) {
        tmp_exprsMatOther <- exprsMatOther[, cells]
      }else{
        tmp_exprsMatOther <- NULL
      }

      tmp_results <- reRunCLUEY(exprsMatRNA= exprsMatRNA[, cells], exprsMatOther=tmp_exprsMatOther, reducedDims = reducedDims[cells,],
                                reference = reference, corMethod = corMethod, kLimit = subK,
                                previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, minCells = minCells,
                                encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)

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
    cluey_df <- cluey_df[match(colnames(exprsMatRNA), cluey_df$cell_id),]

    score <- FisherZInv(mean(FisherZ(unique(cluey_df$correlation))))
    pvalue <- pnorm(mean(FisherZ(unique(cluey_df$correlation))), lower.tail = F)
    return(list(optimal_K = length(unique(cluey_df$annotation)), annotations = cluey_df, score = round(score, digits = 2), pvalue = round(pvalue, digits = 2)))

  }
  message("Returning results...")
  return(current_result)

}


reRunCLUEY <-  function(exprsMatRNA, exprsMatOther=exprsMatOther, reducedDims, reference, corMethod, kLimit, previousResult, propGenes, recursive, minCells,
                        encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs){

  if(ncol(exprsMatRNA) > minCells){

    if(ncol(reducedDims) > encodingDim2){

      if (unimodal) {

        model <- modelGeneVar(exprsMatRNA)
        hvg_genes <- getTopHVGs(model, n=0.025*nrow(exprsMatRNA))
        hvg_data <- exprsMatRNA[hvg_genes, ]

        reducedDims <- getEncodingMultiModal(list(hvg_data), hiddenDims = 10, encodedDim = encodingDim2, epochs=nEpochs, batchSize = 16)
      } else if (!unimodal) {

        model <- modelGeneVar(exprsMatRNA)
        hvg_genes <- getTopHVGs(model, n=0.025*nrow(exprsMatRNA))
        hvg_data <- exprsMatRNA[hvg_genes, ]

        reducedDims <- getEncodingMultiModal(list(hvg_data, exprsMatOther), hiddenDims = c(50, 10), encodedDim = encodingDim2, epochs=nEpochs, batchSize =16)
      }
    }

    dist_matrix <- as.matrix(dist(reducedDims))
    affinity_matrix <- affinityMatrix(dist_matrix, K = 20, sigma = 0.4)

    k_annotations <- list()
    k_clusters <- list()

    for(k in 2:kLimit){
      cluster_res <- spectralClustering(as.matrix(affinity_matrix),  k)
      clusters <- cluster_res$label
      k_clusters[[k-1]] <- clusters
      k_annotations[[k-1]] <- annotateClusters(exprsMat = exprsMatRNA, clusters = clusters,
                                                corMethod = corMethod, reference = reference, propGenes = propGenes)

    }
    names(k_annotations) <- 2:kLimit
    names(k_clusters) <- 2:kLimit
    optimal_annotations <- findOptimalK(k_annotations, k_clusters)

    cluey_df <- data.frame(cell_id = colnames(exprsMatRNA),
                           spectral_cluster = optimal_annotations$clusters,
                           annotation = gsub("[.]", " ", optimal_annotations[[1]]$annotation[optimal_annotations$clusters]),
                           correlation = optimal_annotations[[1]]$correlation[optimal_annotations$clusters])

    cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
    max_cluster <- max(cluey_df$cluster)
    score <- FisherZInv(mean(FisherZ(unique(cluey_df$correlation))))
    pvalue <- pnorm(mean(FisherZ(unique(cluey_df$correlation))), lower.tail = F)

    current_result <- list(optimal_K = max_cluster, score = round(optimal_annotations$score, digits = 2), annotations = cluey_df, pvalue = round(pvalue, digits = 2))
    tmp_current_result <- current_result
    results <- list()

    if(recursive == TRUE){

      if(max_cluster > 1){
        for(i in 1:max_cluster){
          idx <- which(cluey_df$cluster == i)
          cells <- cluey_df$cell_id[idx]

          if(length(cells) > minCells){
            tmp_current_result[["annotations"]] <- tmp_current_result$annotations[tmp_current_result$annotations$cell_id %in% cells,]

            if (!unimodal) {
              tmp_exprsMatOther <- exprsMatOther[, cells]
            }else{
              tmp_exprsMatOther <- NULL
            }

            tmp_results <- reRunCLUEY(exprsMatRNA = exprsMatRNA[, cells], exprsMatOther = tmp_exprsMatOther, reducedDims = reducedDims[cells,],
                                      reference = reference, corMethod = corMethod, kLimit = kLimit,
                                      previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, minCells = minCells,
                                      encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)

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
        cluey_df <- cluey_df[match(colnames(exprsMatRNA), cluey_df$cell_id),]

        score <- FisherZInv(mean(FisherZ(unique(cluey_df$correlation))))
        pvalue <- pnorm(mean(FisherZ(unique(cluey_df$correlation))), lower.tail = F)
        return(list(optimal_K = length(unique(cluey_df$annotation)), annotations = cluey_df, score = round(score, digits = 2), pvalue = round(pvalue, digits = 2)))

      }

    }else{
      message("...")
      return(current_result)
    }
  }
  message("...")
  return(previousResult)
}

generatePseudobulk <- function(exprsMat, clusters){

  cluster_labels <- unique(clusters)
  num_clusters <- length(cluster_labels)
  pseudobulk_profiles <- vector("list", num_clusters)

  for(i in 1:num_clusters){

    profile <- exprsMat[, which(clusters == cluster_labels[i])]

    if(!is.null(ncol(profile))){

      pseudobulk_profiles[[i]] <- rowMeans(exprsMat[, which(clusters == cluster_labels[i])])

    }else{

      pseudobulk_profiles[[i]] <- profile

    }


  }

  names(pseudobulk_profiles) <- cluster_labels
  pseudobulk_profiles
}

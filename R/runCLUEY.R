#' runCLUEY
#'
#' @description
#' Perform cluster estimation of single or multi-modal single-cell data
#'
#' @param exprsMatRNA A numeric matrix of normalised gene expression data where rows are
#' features and columns are cells.
#' @param exprsMatOther A list of additional modalities, such as ATAC-seq or CITE-seq,
#' where rows are the features and columns are cells.
#' @param knowledgeBase knowledgeBase generated from the function `generateKnowledgeBase`.
#' @param corMethod Correlation method used when determining the most similar cell-types
#' for each cluster. Options are "Spearman" or "Pearson". Default is Spearman correlation.
#' @param propGenes Proportion of genes to use when calculating the correlation of the expression
#' profiles between each cluster and knowledgeBase cell-type. Default is 0.25.
#' @param kLimit The maximum number of clusters to estimate. Default is 20.
#' @param subK The maximum number of clusters to estimate when performing iterative clustering. Default is 3.
#' @param minCells The minimum number of cells that a cluster should contain. Default is 20 cells.
#' @param recursive A boolean specifying whether to perform clustering recursively. Default is TRUE.
#' @param encodingDim1 Dimension of the hidden layer for autoencoder at the first level of cluster estimation. Deault is 50.
#' @param encodingDim2 Dimension of the hidden layer for autoencoder at the second level of cluster estimation. Deault is 20.
#' @param hiddenDimsMultimodal Dimension of the hidden layer for multi-modal autoencoder. Ignored if `exprsMatOther` is NULL. Default is 50.
#' @param nEpochs Number of epochs when training the autoencoder.
#'
#' @return A list containing the optimal number of clusters `optimal_k` and
#' final clustering score `score`. Cluster assignments and correlation scores
#' for each cluster can be in `predictions`.
#' @export
#'
runCLUEY <-  function(exprsMatRNA, exprsMatOther=NULL, knowledgeBase, corMethod = "spearman", propGenes=0.25, kLimit=20, subK=3, minCells=20,
                      recursive=TRUE, encodingDim1=50, encodingDim2=10, hiddenDimsMultimodal=50, nEpochs=30){
  unimodal <- TRUE
  stopifnot(minCells >= 2*subK)
  rownames(exprsMatRNA) <- toupper(rownames(exprsMatRNA))

  if (!is.null(exprsMatOther)) {
    unimodal <-  FALSE
    rownames(exprsMatOther) <- toupper(rownames(exprsMatOther))
  }

  if (unimodal) {

    message("Selecting HVGs...\n")
    model <- scran::modelGeneVar(exprsMatRNA)
    hvg_genes <- scran::getTopHVGs(model, n=0.05*nrow(exprsMatRNA))
    hvg_data <- exprsMatRNA[hvg_genes, ]
    message("Reducing dimensions...\n")
    reducedDims <- getEncodingMultiModal(list(hvg_data), hiddenDims = encodingDim1, encodedDim = encodingDim1, epochs=nEpochs, batchSize=32)
  } else if (!unimodal) {
    message("Selecting HVGs...\n")
    model <- scran::modelGeneVar(exprsMatRNA)
    hvg_genes <- scran::getTopHVGs(model, n=0.05*nrow(exprsMatRNA))
    hvg_data <- exprsMatRNA[hvg_genes, ]
    message("Reducing dimensions...\n")
    reducedDims <- getEncodingMultiModal(list(hvg_data, exprsMatOther), hiddenDims = hiddenDimsMultimodal, encodedDim = encodingDim1, epochs=nEpochs, batchSize=32)
  }


  dist_matrix <- as.matrix(distances::distances(reducedDims))
  affinity_matrix <- SNFtool::affinityMatrix(dist_matrix, K = 30, sigma = 0.4)

  k_predictions <- list()
  k_clusters <- list()
  message("Estimating initial k...\n")
  for(k in 2:kLimit){
    cluster_res <- spectralClustering(as.matrix(affinity_matrix),  k)
    clusters <- cluster_res$label
    k_clusters[[length(k_clusters) + 1]] <- clusters
    k_predictions[[length(k_predictions) + 1]] <- annotateClusters(exprsMatRNA, clusters, corMethod = corMethod, knowledgeBase = knowledgeBase, propGenes = propGenes)

  }
  names(k_predictions) <- 2:kLimit
  names(k_clusters) <- 2:kLimit
  optimal_predictions <- findOptimalK(k_predictions, k_clusters) # Pass in dis matrix

  cluey_df <- data.frame(cell_id = colnames(exprsMatRNA),
                         spectral_cluster = optimal_predictions$clusters,
                         annotation = gsub("[.]", " ", optimal_predictions[[1]]$annotation[optimal_predictions$clusters]),
                         correlation = optimal_predictions[[1]]$correlation[optimal_predictions$clusters])

  cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
  cluey_df <- cluey_df[match(colnames(exprsMatRNA), cluey_df$cell_id),]
  max_cluster <- max(cluey_df$cluster)
  score <- DescTools::FisherZInv(mean(DescTools::FisherZ(unique(cluey_df$correlation))))
  pvalue <- pnorm(mean(DescTools::FisherZ(unique(cluey_df$correlation))), lower.tail = F)

  current_result <- list(optimal_K = max_cluster, score = round(optimal_predictions$score, digits = 2), predictions = cluey_df, pvalue = round(pvalue, digits = 2))
  tmp_current_result <- current_result
  results <- list()

  message("Running multiscale clustering...\n")
  for(i in 1:max_cluster){

    idx <- which(cluey_df$cluster == i)
    cells <- cluey_df$cell_id[idx]
    if(length(cells) > minCells){

      tmp_current_result[["predictions"]] <- tmp_current_result$predictions[tmp_current_result$predictions$cell_id %in% cells,]

      if (!unimodal) {
        tmp_exprsMatOther <- exprsMatOther[, cells]
      }else{
        tmp_exprsMatOther <- NULL
      }

      tmp_results <- reRunCLUEY(exprsMatRNA= exprsMatRNA[, cells], exprsMatOther=tmp_exprsMatOther, reducedDims = reducedDims[cells,],
                                knowledgeBase = knowledgeBase, corMethod = corMethod, kLimit = subK,
                                previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, minCells = minCells,
                                encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)

      results[[length(results)+1]] <- tmp_results
    }
  }


  if(length(results) > 0){

    cluey_df <- current_result$predictions
    for(i in 1:length(results)){
      result <- results[[i]]
      if(!is.null(result)){
        cells <- result$predictions$cell_id
        cluey_df <- cluey_df[!(cluey_df$cell_id %in% cells),]
        cluey_df <- rbind(cluey_df, results[[i]]$predictions)
      }

    }
    cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
    cluey_df <- cluey_df[match(colnames(exprsMatRNA), cluey_df$cell_id),]

    score <- DescTools::FisherZInv(mean(DescTools::FisherZ(unique(cluey_df$correlation))))
    pvalue <- pnorm(mean(DescTools::FisherZ(unique(cluey_df$correlation))), lower.tail = F)
    return(list(optimal_K = length(unique(cluey_df$annotation)), predictions = cluey_df, score = round(score, digits = 2), pvalue = round(pvalue, digits = 2)))

  }
  message("Returning results...")
  return(current_result)

}

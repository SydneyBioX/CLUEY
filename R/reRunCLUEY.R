###################################################################
# Function called during recursive clustering stage
###################################################################

reRunCLUEY <-  function(exprsMatRNA, exprsMatOther=exprsMatOther, reducedDims, knowledgeBase, corMethod, kLimit, previousResult, propGenes, recursive, minCells,
                        encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs){

  if(ncol(exprsMatRNA) > minCells){

    if(ncol(reducedDims) > encodingDim2){

      if (unimodal) {

        model <- scran::modelGeneVar(exprsMatRNA)
        hvg_genes <- scran::getTopHVGs(model, n=0.025*nrow(exprsMatRNA))
        hvg_data <- exprsMatRNA[hvg_genes, ]

        reducedDims <- getEncodingMultiModal(list(hvg_data), hiddenDims = 10, encodedDim = encodingDim2, epochs=nEpochs, batchSize = 16)
      } else if (!unimodal) {

        model <- scran::modelGeneVar(exprsMatRNA)
        hvg_genes <- scran::getTopHVGs(model, n=0.025*nrow(exprsMatRNA))
        hvg_data <- exprsMatRNA[hvg_genes, ]

        reducedDims <- getEncodingMultiModal(list(hvg_data, exprsMatOther), hiddenDims = c(50, 10), encodedDim = encodingDim2, epochs=nEpochs, batchSize =16)
      }
    }

    dist_matrix <- as.matrix(dist(reducedDims))
    affinity_matrix <- SNFtool::affinityMatrix(dist_matrix, K = 20, sigma = 0.4)

    k_predictions <- list()
    k_clusters <- list()

    for(k in 2:kLimit){
      cluster_res <- spectralClustering(as.matrix(affinity_matrix),  k)
      clusters <- cluster_res$label
      k_clusters[[k-1]] <- clusters
      k_predictions[[k-1]] <- annotateClusters(exprsMat = exprsMatRNA, clusters = clusters,
                                               corMethod = corMethod, knowledgeBase = knowledgeBase, propGenes = propGenes)

    }
    names(k_predictions) <- 2:kLimit
    names(k_clusters) <- 2:kLimit
    optimal_predictions <- findOptimalK(k_predictions, k_clusters)

    cluey_df <- data.frame(cell_id = colnames(exprsMatRNA),
                           spectral_cluster = optimal_predictions$clusters,
                           annotation = gsub("[.]", " ", optimal_predictions[[1]]$annotation[optimal_predictions$clusters]),
                           correlation = optimal_predictions[[1]]$correlation[optimal_predictions$clusters])

    cluey_df$cluster <- as.integer(factor(cluey_df$annotation))
    max_cluster <- max(cluey_df$cluster)
    score <- DescTools::FisherZInv(mean(DescTools::FisherZ(unique(cluey_df$correlation))))
    pvalue <- pnorm(mean(DescTools::FisherZ(unique(cluey_df$correlation))), lower.tail = F)

    current_result <- list(optimal_K = max_cluster, score = round(optimal_predictions$score, digits = 2), predictions = cluey_df, pvalue = round(pvalue, digits = 2))
    tmp_current_result <- current_result
    results <- list()

    if(recursive == TRUE){

      if(max_cluster > 1){
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

            tmp_results <- reRunCLUEY(exprsMatRNA = exprsMatRNA[, cells], exprsMatOther = tmp_exprsMatOther, reducedDims = reducedDims[cells,],
                                      knowledgeBase = knowledgeBase, corMethod = corMethod, kLimit = kLimit,
                                      previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, minCells = minCells,
                                      encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)

            results[[length(results)+1]] <- tmp_results


          }
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

    }else{
      message("...")
      return(current_result)
    }
  }
  message("...")
  return(previousResult)
}

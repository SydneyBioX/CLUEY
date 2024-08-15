require(tensorflow)
require(luz)
require(SingleCellExperiment)
require(keras)
require(parallel)
require(cluster)
require(stats)
require(DescTools)
require(RANN)
require(SNFtool)
require(scran)
require(Seurat)
require(Cepo)
require(distances)

# Cluey with weighted average correlation coefficients - Fisher's Z method
annotate_clusters <- function(exprsMat, clusters, corMethod, reference, propGenes){
  
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

find_optimal_k <- function(k_annotations, k_clusters, threshold){
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


runCluey <-  function(exprsMatRNA, exprsMatOther=NULL, reference, corMethod = "spearman", propGenes=0.25, kLimit=20, subK=3, minCells=20,
                      recursive=TRUE, threshold=0.5, encodingDim1=50, encodingDim2=10, hiddenDimsMultimodal=50, nEpochs=50){
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
    k_annotations[[length(k_annotations) + 1]] <- annotate_clusters(exprsMatRNA, clusters, corMethod = corMethod, reference = reference, propGenes = propGenes)
    
  }
  names(k_annotations) <- 2:kLimit
  names(k_clusters) <- 2:kLimit
  optimal_annotations <- find_optimal_k(k_annotations, k_clusters, threshold) # Pass in dis matrix
  
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
      
      tmp_results <- reRunCluey(exprsMatRNA= exprsMatRNA[, cells], exprsMatOther=tmp_exprsMatOther, reducedDims = reducedDims[cells,], 
                                reference = reference, corMethod = corMethod, kLimit = subK, 
                                previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, threshold = threshold,
                                minCells = minCells, encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)
      
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


reRunCluey <-  function(exprsMatRNA, exprsMatOther=exprsMatOther, reducedDims, reference, corMethod, kLimit, previousResult, propGenes, recursive, threshold, minCells,
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
      k_annotations[[k-1]] <- annotate_clusters(exprsMat = exprsMatRNA, clusters = clusters,
                                                corMethod = corMethod, reference = reference, propGenes = propGenes)
      
    }
    names(k_annotations) <- 2:kLimit
    names(k_clusters) <- 2:kLimit
    optimal_annotations <- find_optimal_k(k_annotations, k_clusters, threshold) # Pass in dis matrix
    
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
            
            tmp_results <- reRunCluey(exprsMatRNA = exprsMatRNA[, cells], exprsMatOther = tmp_exprsMatOther, reducedDims = reducedDims[cells,], 
                                      reference = reference, corMethod = corMethod, kLimit = kLimit, 
                                      previousResult = tmp_current_result, propGenes = propGenes, recursive = recursive, threshold = threshold,
                                      minCells = minCells, encodingDim2=encodingDim2, unimodal=unimodal, hiddenDimsMultimodal=hiddenDimsMultimodal, nEpochs=nEpochs)
            
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




spectralClustering <- function(affinity, 
                               K = 20, 
                               delta = 1e-5) {
  
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- affinity
  
  
  neff <- K + 1
  
  
  v <- sqrt(d)
  NL <- L/(v %*% t(v))
  
  
  f <- function(x, A = NULL){ 
    as.matrix(A %*% x)
  }
  
  n <- nrow(affinity)
  
  
  NL <- ifelse(NL > delta, NL, 0)
  ###########
  NL <- methods::as(NL, "dgCMatrix")
  
  
  eig <- igraph::arpack(f, extra = NL, sym = TRUE,
                        options = list(which = 'LA', nev = neff,
                                       n = n,
                                       ncv = max(min(c(n,4*neff))),
                                       maxiter=10000000
                        )
  )
  
  
  
  psi <- eig$vectors / (eig$vectors[,1] %*% matrix(1, 1, neff))#right ev
  eigenvals <- eig$values
  
  res <- sort(abs(eigenvals), index.return = TRUE, decreasing = TRUE)
  U <- eig$vectors[, res$ix[seq_len(K)]]
  normalize <- function(x) x/sqrt(sum(x^2))
  
  U <- t(apply(U, 1, normalize))
  eigDiscrete <- .discretisation(U)
  eigDiscrete <- eigDiscrete$discrete
  labels <- apply(eigDiscrete, 1, which.max)
  
  lambda <- eigenvals[-1]/(1 - eigenvals[-1])
  lambda <- rep(1,n) %*% t(lambda)
  lam <- lambda[1,]/lambda[1,1]
  neigen <- K
  eigenvals <- eigenvals[seq_len((neigen + 1))]
  X <- psi[,2:(neigen + 1)]*lambda[, seq_len(neigen)] 
  
  return(list(labels = labels,
              eigen_values = eig$values,
              eigen_vectors = eig$vectors,
              X = X))
}

standardise <- function(x){(x-mean(x))/(sd(x))}

maxmin <- function(x){(x - min(x))/(max(x) - min(x))}

createEncoder <- function(layer_input, output_dim) {
  
  encoder <- layer_input %>%
    layer_dense(units = output_dim) %>% 
    layer_activation_relu() %>% 
    layer_batch_normalization() %>% 
    layer_activation_leaky_relu() 
  # layer_batch_normalization()
  
  return(encoder)
}

createDecoder = function(latent_space, output_dim) {
  
  decoder = latent_space %>%
    layer_dense(units = output_dim) %>% 
    layer_activation_relu() %>%
    layer_batch_normalization() %>%
    layer_activation_leaky_relu() 
  # layer_batch_normalization() %>% 
  # layer_activation_leaky_relu() %>% 
  # layer_batch_normalization()
  
  
  return(decoder)
}

getEncodingMultiModal <- function(datasetList, learningRate = 0.001,
                                  hiddenDims = 50, encodedDim = 10,
                                  epochs = 50, batchSize = 32,
                                  dropoutRate = 0, verbose = 0) {
  
  
  datasetList_t = lapply(datasetList, function(x){
    x <- x[!(rowSums(x) == 0),]
    sds <- colSds(x)
    means <- colMeans(x)
    standardized_matrix <- sweep(x, 2, means, "-")
    standardized_matrix <- sweep(standardized_matrix, 2, sds, "/")
    t(standardized_matrix)
  })
  
  input_dim = lapply(datasetList_t, ncol)
  input_dat = mapply(function(datm, dim) {
    
    list(layer_input(shape = ncol(datm)), dim) 
    
  }, datm = datasetList_t, dim = hiddenDims, SIMPLIFY = F)
  
  datasetList = datasetList_t
  
  # Create encoder for each modality
  each_modal_encoder = lapply(input_dat, function(x){
    
    createEncoder(layer_input = x[[1]], output_dim = x[[2]])
    
  })
  
  # Concatenate encoders
  joint_modals = layer_concatenate(each_modal_encoder)
  
  if(length(datasetList) > 1){
    
    joint_latent_space = joint_modals %>%
      layer_dense(units = encodedDim) %>% 
      layer_activation_leaky_relu() %>% 
      layer_batch_normalization() %>% 
      layer_activation_leaky_relu()
    
  }else{
    
    joint_latent_space = joint_modals
    
  }
  
  decoder_dat = mapply(function(latent_space, dim) {
    
    list(latent_space, dim)
    
  }, latent_space = rep(list(joint_latent_space), times = length(datasetList)), dim = input_dim, SIMPLIFY = F)
  
  each_modal_decoder = lapply(decoder_dat, function(x) {
    
    createDecoder(latent_space = x[[1]], output_dim = x[[2]])
    
  })
  
  autoencoder_multimodal_model <- keras_model(
    inputs = lapply(input_dat,function(x) x[[1]]),
    outputs = each_modal_decoder
  )
  
  ADAM = optimizer_adam(learning_rate = learningRate)
  autoencoder_multimodal_model %>% compile(
    loss = 'mean_squared_error',
    optimizer = ADAM
  )
  
  
  history = autoencoder_multimodal_model %>%
    keras::fit(
      datasetList,
      datasetList,
      epochs = epochs,
      shuffle = TRUE,
      batch_size = batchSize,
      verbose = verbose)
  
  
  
  latent_space_model = keras_model(
    inputs = lapply(input_dat,function(x) x[[1]]),
    outputs = joint_latent_space  # Joint latent space
  )
  
  # Use the new model to predict the joint latent space representation
  joint_latent_representation = latent_space_model %>% predict(datasetList, batch_size = batchSize)
  
  rownames(joint_latent_representation) = rownames(datasetList[[1]])
  
  return(joint_latent_representation)
  
}

.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors <- t(apply(eigenVectors,1,normalize))
  
  n <- nrow(eigenVectors)
  k <- ncol(eigenVectors)
  
  R <- matrix(0,k,k)
  R[,1] <- t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i <- which(x == min(x))
    return(i[1])
  }
  
  c <- matrix(0, n, 1)
  for (j in seq(2, k)) {
    c <- c + abs(eigenVectors %*% matrix(R[,j - 1], k, 1))
    i <- mini(c)
    R[,j] <- t(eigenVectors[i,])
  }
  
  lastObjectiveValue <- 0
  for (i in seq_len(1000)) {
    eigenDiscrete <- .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde <- svd(t(eigenDiscrete) %*% eigenVectors)
    U <- svde[['u']]
    V <- svde[['v']]
    S <- svde[['d']]
    
    NcutValue <- 2 * (n - sum(S))
    if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break
    
    lastObjectiveValue <- NcutValue
    R <- V %*% t(U)
    
  }
  
  return(list(discrete = eigenDiscrete, continuous = eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y <- matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i <- which(x == max(x))
    return(i[1])
  }
  j <- apply(eigenVector,1,maxi)
  Y[cbind(seq_len(nrow(eigenVector)), j)] <- 1
  
  return(Y)
  
}

generate_pseudobulk <- function(exprsMat, clusters){
  
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

generate_reference <- function(reference, celltypes, batch = NULL, method = "ds", minCells = 20){
  if(!(2*ncol(reference) == sum(ncol(reference), length(celltypes)))){
    stop("number of cells and length of celltypes vector do not match. Batch vector must also be the same length if it is being used...")
  }
  
  rownames(reference) <- toupper(rownames(reference))
  celltypes <- make.names(celltypes)
  idx <- rowSums(reference)
  idx <- which(idx == 0)
  if(length(idx) > 0){
    reference <- suppressWarnings(as.matrix(reference[-idx,]))
  }else{
    reference <- suppressWarnings(as.matrix(reference))
  }
  
  if(method == "ds"){
    if(!is.null(batch)){
      message("Identifying differentially stable genes per batch...\n")
      batches <- unique(batch)
      batch_genes <- lapply(batches, function(x){
        idx <- which(batch == x)
        
        n_celltypes <- unique(celltypes[idx])
        if(length(n_celltypes) > 1){
          keep_celltypes <- names(table(celltypes[idx]) >= minCells)
          if(length(keep_celltypes) > 0){
            idx <- which(celltypes %in% keep_celltypes)
            reference_cepo <- suppressWarnings(Cepo(reference[,idx], celltypes[idx], minCells = minCells, exprsPct = 0.05))
            reference_top_genes <- topGenes(reference_cepo, n = nrow(reference_cepo$stats), returnValues = T)
            ref_genes <- lapply(reference_top_genes, function(x){x[x > 0]})
            ref_genes
          }
          
        }
        
      })
      batch_genes <- unlist(batch_genes, recursive = FALSE)
      
      grouped_celltypes <- split(batch_genes, names(batch_genes))
      ref_genes <- lapply(grouped_celltypes, function(x){
        genes <- names(unlist(unname(x)))
        common_genes <- genes[which(duplicated(genes) == TRUE)]
        if(length(common_genes) > 0){
          average_values <- tapply(unlist(unname(x)), names(unlist(unname(x))), mean)
          average_values[order(average_values, decreasing = TRUE)]
        }else{
          unlist(unname(x))
        }
        
      })
    }else{
      message("Identifying differentially stable genes...\n")
      if(sum(table(celltypes) < minCells) > 0){
        stop(paste0("Less than ", minCells, " cells for at least one cell-type, please filter out or decrease mincells parameter..."))
      }
      
      reference_cepo <- suppressWarnings(Cepo(reference, celltypes, minCells = minCells, exprsPct = 0.05))
      reference_top_genes <- topGenes(reference_cepo, n = nrow(reference_cepo$stats), returnValues = T)
      ref_genes <- lapply(reference_top_genes, function(x){x[x > 0]})
    }
    
    
  }
  
  if(method == "de"){
    message("Identifying differentially expressed genes...\n")
    seurat_ref <- suppressWarnings(CreateSeuratObject(counts = reference, assay = "RNA", meta.data = data.frame(celltypes = celltypes)))
    seurat_ref <- NormalizeData(seurat_ref)
    Idents(object = seurat_ref) <- celltypes
    tmp <- FindAllMarkers(seurat_ref, test.use="t")
    ref_genes <- lapply(unique(tmp$cluster), function(celltype){
      ct_tmp <- tmp[tmp$cluster %in% celltype,]
      ct_tmp <- ct_tmp[order(ct_tmp$p_val_adj, decreasing = FALSE),]
      ct_tmp[ct_tmp$p_val_adj < 0.05,]$gene
    })
    names(ref_genes) <- make.names(levels(seurat_ref))
  }
  
  reference_pseudobulk <- mapply(function(genes, celltype){
    if(method == "ds"){
      genes <- names(genes)
    }else{
      genes <- intersect(genes, rownames(reference))
    }
    tmp <- reference[genes, which(celltypes == celltype)]
    tmp <- rowMeans(tmp)
    tmp
    
  },ref_genes, names(ref_genes))
  
  
  
  return(reference_pseudobulk)
}


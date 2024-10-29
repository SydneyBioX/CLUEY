###################################################################
# Function to extract muli-modal latent space for clustering
###################################################################

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

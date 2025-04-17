###################################################################
# Function to extract muli-modal latent space for clustering
###################################################################
#' @importFrom magrittr %>%

get_encoding <- function(dataset_list, learning_rate = 0.001,
                         hidden_dims = 50, encoded_dim = 10,
                         epochs = 50, batch_size = 32,
                         dropoutRate = 0, verbose = 0) {


  dataset_list_t = lapply(dataset_list, function(x){
    x <- x[,!(colSums(x) == 0)]
    sds <- colSds(x)
    means <- colMeans(x)
    standardized_matrix <- sweep(x, 2, means, "-")
    standardized_matrix <- sweep(standardized_matrix, 2, sds, "/")
    standardized_matrix
  })

  input_dim <-  lapply(dataset_list_t, ncol)
  input_dat <-  mapply(function(datm, dim) {

    list(keras::layer_input(shape = ncol(datm)))

  }, datm = dataset_list_t, SIMPLIFY = F)

  dataset_list = dataset_list_t

  # Create encoder for each modality
  each_modal_encoder = lapply(input_dat, function(x){

    create_encoder(layer_input = x[[1]], hidden_dim=hidden_dims, output_dim=encoded_dim)
  })

  # Concatenate encoders
  joint_modals = keras::layer_concatenate(each_modal_encoder)

  if(length(dataset_list) > 1){

    joint_latent_space = joint_modals %>%
      keras::layer_dense(units = encoded_dim) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_activation_leaky_relu()

  }else{

    joint_latent_space = joint_modals

  }

  decoder_dat = mapply(function(latent_space, dim) {

    list(latent_space, dim)

  }, latent_space = rep(list(joint_latent_space), times = length(dataset_list)), dim = input_dim, SIMPLIFY = F)

  each_modal_decoder = lapply(decoder_dat, function(x) {

    create_decoder(latent_space = x[[1]], output_dim = x[[2]])

  })

  autoencoder_multimodal_model <- keras::keras_model(
    inputs = lapply(input_dat, function(x) x[[1]]),
    outputs = each_modal_decoder
  )

  ADAM = keras::optimizer_adam(learning_rate = learning_rate)
  autoencoder_multimodal_model %>% keras::compile(
    loss = 'mean_squared_error',
    optimizer = ADAM
  )


  history = autoencoder_multimodal_model %>%
    keras::fit(
      dataset_list,
      dataset_list,
      epochs = epochs,
      shuffle = TRUE,
      batch_size = batch_size,
      verbose = verbose)

  latent_space_model = keras::keras_model(
    inputs = lapply(input_dat, function(x) x[[1]]),
    outputs = joint_latent_space  # Joint latent space
  )

  # Use the new model to predict the joint latent space representation
  joint_latent_representation = latent_space_model %>% predict(dataset_list, batch_size = batch_size)

  rownames(joint_latent_representation) = rownames(dataset_list[[1]])

  return(joint_latent_representation)

}

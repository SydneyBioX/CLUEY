###################################################################
# Function for creating decoder
###################################################################
#' @importFrom magrittr %>%
createDecoder = function(latent_space, output_dim) {

  decoder = latent_space %>%
    keras::layer_dense(units = output_dim) %>%
    keras::layer_activation_relu() %>%
    keras::layer_batch_normalization() %>%
    keras::layer_activation_leaky_relu()

  return(decoder)
}

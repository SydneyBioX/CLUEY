###################################################################
# Function for creating decoder
###################################################################

createDecoder = function(latent_space, output_dim) {

  decoder = latent_space %>%
    layer_dense(units = output_dim) %>%
    layer_activation_relu() %>%
    layer_batch_normalization() %>%
    layer_activation_leaky_relu()

  return(decoder)
}

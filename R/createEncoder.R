###################################################################
# Function for creating encoder
###################################################################
createEncoder <- function(layer_input, output_dim) {

  encoder <- layer_input %>%
    layer_dense(units = output_dim) %>%
    layer_activation_relu() %>%
    layer_batch_normalization() %>%
    layer_activation_leaky_relu()

  return(encoder)
}

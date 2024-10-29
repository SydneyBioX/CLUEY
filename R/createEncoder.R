###################################################################
# Function for creating encoder
###################################################################
#' @importFrom magrittr %>%
createEncoder <- function(layer_input, output_dim) {

  encoder <- layer_input %>%
    keras::layer_dense(units = output_dim) %>%
    keras::layer_activation_relu() %>%
    keras::layer_batch_normalization() %>%
    keras::layer_activation_leaky_relu()

  return(encoder)
}


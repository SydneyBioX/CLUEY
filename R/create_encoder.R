###################################################################
# Function for creating encoder
###################################################################
#' @importFrom magrittr %>%

create_encoder <- function(layer_input, hidden_dim, output_dim) {

  encoder <- layer_input %>%
    keras::layer_dense(units = hidden_dim) %>%
    keras::layer_activation_leaky_relu() %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = hidden_dim/2) %>%
    keras::layer_activation_leaky_relu() %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = output_dim) %>%
    keras::layer_activation_leaky_relu()

  return(encoder)
}

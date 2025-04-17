###################################################################
# Function for creating decoder
###################################################################
#' @importFrom magrittr %>%
#'
create_decoder = function(latent_space, output_dim) {

    decoder = latent_space %>%
    keras::layer_dense(units = output_dim) %>%
    keras::layer_activation_relu()

  return(decoder)
}

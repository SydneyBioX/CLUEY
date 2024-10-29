###################################################################
# Function for standardising expression matrix
###################################################################
standardise <- function(x){(x-mean(x))/(sd(x))}

###################################################################
# Function for .discretisation of Eigen Vector Data
###################################################################
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

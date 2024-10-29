###################################################################
# Function for discretisation
###################################################################
.discretisation <- function(eigenVectors) {

  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors <- t(apply(eigenVectors,1,normalize))

  n <- nrow(eigenVectors)
  k <- ncol(eigenVectors)

  R <- matrix(0,k,k)
  R[,1] <- t(eigenVectors[round(n/2),])

  mini <- function(x) {
    i <- which(x == min(x))
    return(i[1])
  }

  c <- matrix(0, n, 1)
  for (j in seq(2, k)) {
    c <- c + abs(eigenVectors %*% matrix(R[,j - 1], k, 1))
    i <- mini(c)
    R[,j] <- t(eigenVectors[i,])
  }

  lastObjectiveValue <- 0
  for (i in seq_len(1000)) {
    eigenDiscrete <- .discretisationEigenVectorData(eigenVectors %*% R)

    svde <- svd(t(eigenDiscrete) %*% eigenVectors)
    U <- svde[['u']]
    V <- svde[['v']]
    S <- svde[['d']]

    NcutValue <- 2 * (n - sum(S))
    if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue <- NcutValue
    R <- V %*% t(U)

  }

  return(list(discrete = eigenDiscrete, continuous = eigenVectors))
}

###################################################################
# Function for performing spectral clustering
###################################################################
spectralClustering <- function(affinity,
                               K = 20,
                               delta = 1e-5) {

  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- affinity


  neff <- K + 1


  v <- sqrt(d)
  NL <- L/(v %*% t(v))


  f <- function(x, A = NULL){
    as.matrix(A %*% x)
  }

  n <- nrow(affinity)


  NL <- ifelse(NL > delta, NL, 0)
  ###########
  NL <- methods::as(NL, "dgCMatrix")


  eig <- igraph::arpack(f, extra = NL, sym = TRUE,
                        options = list(which = 'LA', nev = neff,
                                       n = n,
                                       ncv = max(min(c(n,4*neff))),
                                       maxiter=10000000
                        )
  )



  psi <- eig$vectors / (eig$vectors[,1] %*% matrix(1, 1, neff))#right ev
  eigenvals <- eig$values

  res <- sort(abs(eigenvals), index.return = TRUE, decreasing = TRUE)
  U <- eig$vectors[, res$ix[seq_len(K)]]
  normalize <- function(x) x/sqrt(sum(x^2))

  U <- t(apply(U, 1, normalize))
  eigDiscrete <- .discretisation(U)
  eigDiscrete <- eigDiscrete$discrete
  labels <- apply(eigDiscrete, 1, which.max)

  lambda <- eigenvals[-1]/(1 - eigenvals[-1])
  lambda <- rep(1,n) %*% t(lambda)
  lam <- lambda[1,]/lambda[1,1]
  neigen <- K
  eigenvals <- eigenvals[seq_len((neigen + 1))]
  X <- psi[,2:(neigen + 1)]*lambda[, seq_len(neigen)]

  return(list(labels = labels,
              eigen_values = eig$values,
              eigen_vectors = eig$vectors,
              X = X))
}

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

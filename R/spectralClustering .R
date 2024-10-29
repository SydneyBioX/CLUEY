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

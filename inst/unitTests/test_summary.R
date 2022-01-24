test_summary <- function() {

  ranks  <- 5:8
  kvec   <- 1:10
  nranks <- length(ranks) 
  nk     <- length(kvec)
  N      <- nranks*nk
  v1     <- rep(ranks, each=nk)
  v2     <- rep(kvec, times=nranks)
  v3     <- 100 + 1/(1:N)
  v4     <- 100 - 1/(1:N)
  x      <- cbind(v1, v2, 999, v3, v4)

  obj    <- getSummary(x, 300, NR=96)

  checkEquals(obj$rank, 5)
  
}

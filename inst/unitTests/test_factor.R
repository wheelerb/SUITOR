test_factor <- function() {

  set.seed(10, kind="Mersenne-Twister", normal.kind="Inversion")

  x    <- matrix(1:6, nrow=2, ncol=3, byrow=FALSE)
  m1   <- matrix(c(1.285714, 3, 4.714286, 1.714286, 4, 6.285714), nrow=2, byrow=TRUE)
  W1   <- matrix(c(2.093777, 2.791703), nrow=2)
  H1   <- matrix(c(0.6140645, 1.432817, 2.25157), nrow=1)
  m2   <- matrix(1:6, nrow=2, byrow=FALSE)
  W2   <- matrix(c(1.637889, 1.509121, 1.452761, 3.624534), nrow=2, byrow=TRUE)
  H2   <- matrix(c(0.1619284, 1.2918984, 2.4218684, 0.4868920, 0.5857803, 0.6846686), nrow=2, byrow=TRUE)

  res1 <- SUITOR:::my_nmf_C(x, 1)
  res2 <- SUITOR:::my_nmf_C(x, 2)

  checkEqualsNumeric(m1, res1$mat, tolerance=1.0e-4)
  checkEqualsNumeric(W1, res1$W,   tolerance=1.0e-4)
  checkEqualsNumeric(H1, res1$H,   tolerance=1.0e-4)
  checkEqualsNumeric(m2, res2$mat, tolerance=1.0e-4)
  checkEqualsNumeric(W2, res2$W,   tolerance=1.0e-4)
  checkEqualsNumeric(H2, res2$H,   tolerance=1.0e-4)

}

cat(paste("\n Starting raoI tests", date(), "\n"))

test_that("checking numerical correction of the raoI function", {
  ngenes <- 7
  ## g1 = WT, d(g1, g2) = d(g1, g3) = 1, d(g2, g3) = 2
  g1 <- g2 <- g3 <- rep(0, ngenes)
  g2[1] <- 1
  g3[2] <- 1
  
  Genotypes1 <- t(t(g1))
  Genotypes2 <- cbind(Genotypes1, t(t(g2)))
  Genotypes3 <- cbind(Genotypes2, t(t(g3)))
  
  Dmat1 <- 0 ## If you pass dist() a vector instead of a matrix, it will 
             ## calculate the distance between components of the vector
             ## instead of the vector with itself (which of course is zero),
             ## as opossed to distance between rows when matrix is given.
  Dmat2 <- as.matrix(dist(t(Genotypes2), method = "manhattan"))
  Dmat3 <- as.matrix(dist(t(Genotypes3), method = "manhattan"))
  
  totalpop <- 500
  
  pops1 <- 500
  pops2 <- c(250, 250)
  pops31 <- c(250, 150, 100)
  pops32 <- c(150, 100, 250)
  
  raoI1 <- OncoSimulR:::raoI(pops1, Dmat1)
  raoI2 <- OncoSimulR:::raoI(pops2, Dmat2)
  raoI31 <- OncoSimulR:::raoI(pops31, Dmat3)
  raoI32 <- OncoSimulR:::raoI(pops32, Dmat3)
  
  raoI1check <- 0
  raoI2check <- 0.5 * 1 * 0.5 + 0.5 * 1 * 0.5
  raoI31check <- (0.5 * (1 * 0.3 + 1 * 0.2)
                  + 0.3 * (1 * 0.5 + 2 * 0.2)
                  + 0.2 * (1 * 0.5 + 2 * 0.3))
  raoI32check <- (0.3 * (1 * 0.2 + 1 * 0.5)
                  + 0.2 * (1 * 0.3 + 2 * 0.5)
                  + 0.5 * (1 * 0.3 + 2 * 0.2))

  expect_equal(raoI1, raoI1check)
  expect_equal(raoI2, raoI2check)
  expect_equal(raoI31, raoI31check)
  expect_equal(raoI32, raoI32check)
})

test_that("checking numerical correction of the raoI function, more", {
  ngenes <- 7
  ## g1 = WT, d(g1, g2) = 2, d(g1, g3) = 7, d(g2, g3) = 5
  g1 <- g2 <- rep(0, ngenes)
  g2[1] <- g2[2] <- 1
  g3 <- rep(1, ngenes)
  
  Genotypes1 <- t(t(g1))
  Genotypes2 <- cbind(Genotypes1, t(t(g2)))
  Genotypes3 <- cbind(Genotypes2, t(t(g3)))
  
  Dmat1 <- as.matrix(dist(t(Genotypes1), method = "manhattan"))
  Dmat2 <- as.matrix(dist(t(Genotypes2), method = "manhattan"))
  Dmat3 <- as.matrix(dist(t(Genotypes3), method = "manhattan"))
  
  pops1 <- 500
  pops2 <- c(250, 250)
  pops31 <- c(250, 150, 100)
  pops32 <- c(150, 100, 250)
  
  raoI1 <- OncoSimulR:::raoI(pops1, Dmat1)
  raoI2 <- OncoSimulR:::raoI(pops2, Dmat2)
  raoI31 <- OncoSimulR:::raoI(pops31, Dmat3)
  raoI32 <- OncoSimulR:::raoI(pops32, Dmat3)
  
  raoI1check <- 0
  raoI2check <- 0.5 * 2 * 0.5 + 0.5 * 2 * 0.5
  raoI31check <- (0.5 * (2 * 0.3 + 7 * 0.2)
                  + 0.3 * (2 * 0.5 + 5 * 0.2)
                  + 0.2 * (7 * 0.5 + 5 * 0.3))
  raoI32check <- (0.3 * (2 * 0.2 + 7 * 0.5)
                  + 0.2 * (2 * 0.3 + 5 * 0.5)
                  + 0.5 * (7 * 0.3 + 5 * 0.2))
  
  expect_equal(raoI1, raoI1check)
  expect_equal(raoI2, raoI2check)
  expect_equal(raoI31, raoI31check)
  expect_equal(raoI32, raoI32check)
})

test_that("checking numerical correction of the raoI function, random", {
  ngenes <- 10
  ngenotypes <- 20
  Genotypes <- matrix(0, nrow = ngenes, ncol = ngenotypes)
  ## First genotype will be WT.
  ## The last genotype will be that with all genes mutated.
  ## The rest will be random, only checking none are duplicated
  while(any(duplicated(t(Genotypes)))) {
    Genotypes[, 2:(ngenotypes-1)] <- matrix(
      sample(0:1, ngenes * (ngenotypes - 2), replace = TRUE),
      nrow = ngenes, ncol = ngenotypes - 2)
    Genotypes[, ngenotypes] <- rep(1, ngenes)
  }
  
  ## Alternative way to calculate Dmat. This makes more calculations
  ## than necessary blablabla, this is just a test. 
  Dmatalt <- apply(Genotypes, 2, function(a) 
    apply(Genotypes, 2, function(b) sum(xor(a, b))))
  Dmat <- as.matrix(dist(t(Genotypes), method = "manhattan"))

  freqs <- sample(0:200, ngenotypes)
  totalpop <- sum(freqs)
  freqsn <- freqs/totalpop
  
  raoI <- OncoSimulR:::raoI(freqs, Dmat)
  raoIcheck <- OncoSimulR:::raoI(freqs, Dmatalt)
  raoIcheck2 <- sum(unlist(
    lapply(1:ngenotypes, function(n)
      lapply(1:ngenotypes, function(m) 
        freqsn[n] * sum(xor(Genotypes[, n], Genotypes[, m])) * freqsn[m]))))
  
  expect_equivalent(Dmat, Dmatalt)
  expect_equal(raoI, raoIcheck)
  expect_equal(raoI, raoIcheck2)
})

test_that("checking numerical correction of the raoI function, errors", {
  ngenes <- 10
  ngenotypes <- 20
  Genotypes <- matrix(0, nrow = ngenes, ncol = ngenotypes)
  ## First genotype will be WT.
  ## The last genotype will be that with all genes mutated.
  ## The rest will be random, only checking none are duplicated
  while(any(duplicated(t(Genotypes)))) {
    Genotypes[, 2:(ngenotypes-1)] <- matrix(
      sample(0:1, ngenes * (ngenotypes - 2), replace = TRUE),
      nrow = ngenes, ncol = ngenotypes - 2)
    Genotypes[, ngenotypes] <- rep(1, ngenes)
  }
  Dmat <- as.matrix(dist(t(Genotypes), method = "manhattan"))
  freqs <- sample(0:200, ngenotypes)

  expect_error(OncoSimulR:::raoI(freqs, "string"), 
               "Unknown argument for dissimilarity matrix")
  expect_error(OncoSimulR:::raoI(freqs, NULL), 
               "Unknown argument for dissimilarity matrix")
  expect_error(OncoSimulR:::raoI(freqs[-1], Dmat), 
               "Dimensions of the dissimilarity matrix are incorrect")
  expect_error(OncoSimulR:::raoI(freqs, Dmat[-1, ]), 
               "Dimensions of the dissimilarity matrix are incorrect")

})

test_that("We obtain same accessible genotypes with different functions", {
    ## More likely to catch bugs if many iters, rather than large matrices
    niter <- 100
    for(i in 1:niter) {
        ## cat("\n i   iteration fast accessible comp ", i)
        for(ng in 2:6) {
            rtmp <- rfitness(ng)
            a1 <- OncoSimulR:::faster_accessible_genotypes_R(rtmp, 0)
            ajm <- OncoSimulR:::genot_to_adj_mat(rtmp)
            a2 <- colnames(OncoSimulR:::filter_inaccessible(ajm, 0))
            a3 <- OncoSimulR:::wrap_accessibleGenotypes(rtmp, 0)
            stopifnot(identical(as.integer(a1), a3))
            stopifnot(identical(as.integer(a2), a3))
            stopifnot(all(a2 ==  a3))
            stopifnot(identical(a2, as.character(a3)))
        }
    } 
})

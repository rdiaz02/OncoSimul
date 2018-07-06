inittime <- Sys.time()
cat(paste("\n Starting accessible_genotypes at", date(), "\n"))
test_that("We obtain same accessible genotypes with different functions plus genot_to_adjm_mat", {
    ## More likely to catch bugs if many iters, rather than large matrices
    niter <- 100
    for(i in 1:niter) {
        ## cat("\n i   iteration fast accessible comp ", i)
        for(ng in 2:6) {
            rtmp <- rfitness(ng)

            ajm <- OncoSimulR:::genot_to_adj_mat(rtmp)
            ajmr <- OncoSimulR:::genot_to_adj_mat_R(rtmp)
            stopifnot(all.equal(ajm, ajmr))
            
            a1 <- OncoSimulR:::faster_accessible_genotypes_R(rtmp, 0)
            a2 <- colnames(OncoSimulR:::filter_inaccessible(ajm, 0))
            a3 <- OncoSimulR:::wrap_accessibleGenotypes(rtmp, 0)
            a4 <- OncoSimulR:::wrap_accessibleGenotypes_former(rtmp, 0)

            stopifnot(identical(as.integer(a1), a3))
            stopifnot(identical(as.integer(a2), a3))
            stopifnot(all(a3 ==  a4))

        }
    } 
})


test_that("The indices of accessible genotypes are correct", {
    ## Make sure we do not assume matrix is ordered

    mf1 <- rbind(c(0, 0, 1),
                 c(1, 0, 4),
                 c(0, 1, .2),
                 c(1, 1, 3))
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf1, 0),
                 c(1, 2))
    mf2 <- rbind(c(0, 0, 1),
                 c(1, 0, 4),
                 c(0, 1, .2),
                 c(1, 1, 5))
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf2, 0),
                 c(1, 2, 4))


    mf1 <- rbind(c(0, 0, 1),
                 c(0, 1, .2),
                 c(1, 0, 4),
                 c(1, 1, 3))
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf1, 0),
                 c(1, 3))
    
    mf2 <- rbind(
        c(0, 1, .2),
        c(0, 0, 1),
        c(1, 0, 4),
        c(1, 1, 5)
    )
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf2, 0),
                 c(2, 3, 4))

    mf2 <- rbind(
        c(0, 1, .2),
        c(0, 0, 1),
        c(1, 1, 5),
        c(1, 0, 4)
    )
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf2, 0),
                 c(2, 3, 4))

    mf2 <- rbind(
        c(0, 1, .2),
        c(0, 0, 1),
        c(1, 1, 5),
        c(1, 0, .4)
    )
    expect_equal(OncoSimulR:::wrap_accessibleGenotypes(mf2, 0),
                 c(2))
    
})


test_that("The indices of adjancey matrices are correct", {
    ## Make sure we do not assume matrix is ordered

    trueam <- matrix(NA, nrow = 4, ncol = 4)
    trueam[1, 2] <- 3
    trueam[1, 3] <- -0.8
    trueam[2, 4] <- -1
    trueam[3, 4] <- 2.8    
    colnames(trueam) <- rownames(trueam) <- 1:4

    mf1 <- rbind(c(0, 0, 1),
                 c(1, 0, 4),
                 c(0, 1, .2),
                 c(1, 1, 3))
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf1),
                          OncoSimulR:::genot_to_adj_mat_R(mf1)))
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf1),
                          trueam[1:4, 1:4]))

    ## these two make use of the fact that we are forced to reorder some
    mf2 <- mf1[c(2, 1, 4, 3), ]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          OncoSimulR:::genot_to_adj_mat_R(mf2)))
    trueam2 <- trueam
    colnames(trueam2) <- rownames(trueam2) <- colnames(trueam)[c(2, 1, 4, 3)]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          trueam2))

    ## ## this is potentially confusing. See below for much cleaner
    ## ## which compares matrices ordered by their names
    ## mf2 <- mf1[c(3, 1, 4, 2), ]
    ## expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
    ##                       OncoSimulR:::genot_to_adj_mat_R(mf2)))
    ## trueam2 <- trueam[c(1, 3, 2, 4), c(1, 3, 2, 4)]
    ## colnames(trueam2) <- rownames(trueam2) <- c(2, 1, 4, 3)
    ## expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
    ##                       trueam2))


    

    ## the next are cleaner:  I compare the matrices ordered by the new names
    ii <- c(2, 1, 4, 3)
    mf2 <- mf1[ii, ]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          OncoSimulR:::genot_to_adj_mat_R(mf2)))
    trueam2 <- trueam[ii, ii]
    colnames(trueam2) <- rownames(trueam2) <- 1:nrow(trueam2)
    ogammf2 <- OncoSimulR:::genot_to_adj_mat(mf2)
    ogammf2 <- ogammf2[order(colnames(ogammf2)), order(colnames(ogammf2))]
    expect_true(all.equal(ogammf2, trueam2))
    expect_true(sum(ogammf2 == trueam2, na.rm = TRUE) == 4)
    expect_false(all(colnames(OncoSimulR:::genot_to_adj_mat(mf2)) == colnames(trueam)))
    expect_false(sum(ogammf2 == trueam, na.rm = TRUE) == 4) ## because 2 and 3 are flipped

    
    

    ii <- c(3, 1, 2, 4)
    mf2 <- mf1[ii, ]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          OncoSimulR:::genot_to_adj_mat_R(mf2)))
    trueam2 <- trueam[ii, ii]
    colnames(trueam2) <- rownames(trueam2) <- 1:nrow(trueam2)
    ogammf2 <- OncoSimulR:::genot_to_adj_mat(mf2)
    ogammf2 <- ogammf2[order(colnames(ogammf2)), order(colnames(ogammf2))]
    expect_true(all.equal(ogammf2, trueam2))
    expect_true(sum(ogammf2 == trueam2, na.rm = TRUE) == 4)
    expect_false(all(colnames(OncoSimulR:::genot_to_adj_mat(mf2)) == colnames(trueam)))
    expect_false(sum(ogammf2 == trueam, na.rm = TRUE) == 4) ## because 2 and 3 are flipped
    
  


    ii <- c(1, 3, 2, 4)
    mf2 <- mf1[ii, ]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          OncoSimulR:::genot_to_adj_mat_R(mf2)))
    trueam2 <- trueam[ii, ii]
    colnames(trueam2) <- rownames(trueam2) <- 1:nrow(trueam2)
    ogammf2 <- OncoSimulR:::genot_to_adj_mat(mf2)
    ogammf2 <- ogammf2[order(colnames(ogammf2)), order(colnames(ogammf2))]
    expect_true(all.equal(ogammf2, trueam2))
    expect_true(sum(ogammf2 == trueam2, na.rm = TRUE) == 4)
    ## note this
    expect_true(all(colnames(OncoSimulR:::genot_to_adj_mat(mf2)) == colnames(trueam)))
    expect_false(sum(ogammf2 == trueam, na.rm = TRUE) == 4) ## because 2 and 3 are flipped
    


    ii <- c(4, 3, 1, 2)
    mf2 <- mf1[ii, ]
    expect_true(all.equal(OncoSimulR:::genot_to_adj_mat(mf2),
                          OncoSimulR:::genot_to_adj_mat_R(mf2)))
    trueam2 <- trueam[ii, ii]
    colnames(trueam2) <- rownames(trueam2) <- 1:nrow(trueam2)
    ogammf2 <- OncoSimulR:::genot_to_adj_mat(mf2)
    ogammf2 <- ogammf2[order(colnames(ogammf2)), order(colnames(ogammf2))]
    expect_true(all.equal(ogammf2, trueam2))
    expect_true(sum(ogammf2 == trueam2, na.rm = TRUE) == 4)
    expect_false(all(colnames(OncoSimulR:::genot_to_adj_mat(mf2)) == colnames(trueam)))
    expect_false(sum(ogammf2 == trueam, na.rm = TRUE) == 4) ## because 2 and 3 are flipped
   
})


## For peaks, see code in test.plotFitnessLandscape.R
cat(paste("\n Ending accessible_genotypes at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

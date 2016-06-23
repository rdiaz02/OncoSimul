test_that("Exercise plotting", {
    r1 <- rfitness(4)
    expect_silent(plot(r1))
    
    ## Specify fitness in a matrix, and plot it
    m5 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 2, 3, 5.5))
    expect_silent(plotFitnessLandscape(m5))
    
    ## Specify fitness with allFitnessEffects, and plot it
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE)))
    
    ## same as
    expect_silent(plotFitnessLandscape(evalAllGenotypes(fe, order = FALSE)))
    
})



test_that("internal peak valley functions", {
    
    x <- matrix(NA, 14, 14)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    
    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)

    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    
    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)



    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 7, 8, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)



    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1 ## hummm.. 3 and 4 should be a peak?Nope, from 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- 1
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(3, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)



    x <- matrix(NA, 5, 5)
    x[1, 3] <- -2
    x[2, 3] <- 4
    x[3, 4] <- 0
    x[4, 5] <- 6

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 5), pv$peak)
    expect_equal(c(2), pv$valley)


    x <- matrix(NA, 5, 5)
    x[1, 3] <- -2
    x[2, 3] <- -5
    x[3, 4] <- 0
    x[4, 5] <- 6

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 2, 5), pv$peak)
    expect_equal(c(3, 4), pv$valley)

    
})

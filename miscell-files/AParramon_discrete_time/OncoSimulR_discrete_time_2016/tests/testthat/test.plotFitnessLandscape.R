test_that("Exercise plotting and dealing with different matrix input", {
    r1 <- rfitness(4)
    expect_silent(plot(r1))
    expect_silent(plot(r1, log = TRUE))
    expect_silent(plot(r1, log = TRUE, use_ggrepel = TRUE))
    expect_silent(plot(r1, log = TRUE, show_labels = FALSE))
    
    
    ## Specify fitness in a matrix, and plot it
    m5 <- cbind(A = c(0, 1, 0, 1), B = c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_silent(plotFitnessLandscape(m5))

    m6 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m6),
                   "Setting/resetting gene names because", fixed = TRUE)

    m7 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m7),
                   "Setting/resetting gene names because", fixed = TRUE)

    m8 <- cbind(A = c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m8),
                   "Setting/resetting gene names because", fixed = TRUE)

    
    ## Specify fitness with allFitnessEffects, and plot it
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE)))

    ## same as
    expect_silent(plotFitnessLandscape(evalAllGenotypes(fe, order = FALSE)))
    ## more ggrepel
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE), use_ggrepel = TRUE))
})


test_that("to_FitnessMatrix stops as it should", {
    x1 <- data.frame(a = 1:2, b = 1:2)
    expect_error(OncoSimulR:::to_Fitness_Matrix(x1, 2000),
                 "We cannot guess what you are passing",
                 fixed = TRUE)
    x2 <- list(a = 12, b = 13)
    expect_error(OncoSimulR:::to_Fitness_Matrix(x2, 2000),
                 "We cannot guess what you are passing",
                 fixed = TRUE)
    ## This is done above
    ## g <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
    ## s1 <- c(1, 1.4, 1.2, 1.5)
    ## expect_error(OncoSimulR:::to_Fitness_Matrix(cbind(g, s1), 2000),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
    ## expect_message(plotFitnessLandscape(cbind(g, s1)),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
    ## expect_message(plotFitnessLandscape(cbind(g, A = c(1, 2))),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
})



test_that("to_FitnessMatrix can deal with df", {
    m4 <- data.frame(G = c("A, B", "A", "WT", "B"),
                     Fitness = c(3, 2, 1, 4))
    expect_message(OncoSimulR:::to_Fitness_Matrix(m4, 2000),
                   "Column names of object", fixed = TRUE)
    m5 <- data.frame(G = c("A, B", "B"),
                     Fitness = c(3, 2))
    expect_message(OncoSimulR:::to_Fitness_Matrix(m5, 2000),
                   "Column names of object", fixed = TRUE)
    x1 <- data.frame(a = c("A, B"), Fitness = 2)
    expect_message(OncoSimulR:::to_Fitness_Matrix(x1, 2000),
                   "Column names of object", fixed = TRUE)
    x2 <- data.frame(a = c("A, B", "B"), Fitness = c(2, 3))
    expect_message(OncoSimulR:::to_Fitness_Matrix(x2, 2000),
                   "Column names of object", fixed = TRUE)
    x3 <- data.frame(a = c("A, B", "C"), Fitness = c(2, 3))
    expect_message(OncoSimulR:::to_Fitness_Matrix(x3, 2000),
                   "Column names of object", fixed = TRUE)
    ## Now, the user code
    expect_message(plotFitnessLandscape(x1))
    expect_message(plotFitnessLandscape(x2))
    expect_message(plotFitnessLandscape(x3))
    expect_message(plotFitnessLandscape(m5))
    expect_message(plotFitnessLandscape(m4))
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

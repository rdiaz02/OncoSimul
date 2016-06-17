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

test_that("mutationPropGrowth warning with Bozic, sample", {
    ## Bozic will crash on purpose if any s >= 1. Make sure
    ## that does not happen here

    RNGkind("Mersenne-Twister")
    set.seed(13)
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = runif(5, 0.01, 0.06),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    
    set.seed(10)
    expect_warning(oiI1 <- oncoSimulSample(4,
                                           oi,
                                           sampleEvery = 0.03,
                                           onlyCancer = FALSE,
                                           detectionProb = NA,
                                           model = "Bozic",
                                           mutationPropGrowth = TRUE,
                                           seed = NULL),
                   "Using fitness Bozic (bozic1) with mutationPropGrowth = TRUE;",
                   fixed = TRUE)
    
})

set.seed(NULL)

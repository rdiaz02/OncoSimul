cat(paste("\n Starting warning-mutPropGrowth tests", date(), "\n"))

## RNGkind("Mersenne-Twister")

test_that("mutationPropGrowth warning with Bozic, indiv", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = runif(5, 0.01, 0.06),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2"),
               drvNames = c("d1", "d2", "f1", "f2", "f3"))
              expect_warning(oiI1 <- oncoSimulIndiv(oi,
                                                    sampleEvery = 0.03,
                                                    keepEvery = 2,
                                                    model = "Bozic",
                                                    mutationPropGrowth = TRUE),
                             "Using fitness Bozic (bozic1) with mutationPropGrowth = TRUE;",
                             fixed = TRUE)                                                     
})


## With seed 5207947 we get a Recoverable exception ti set to
## DBL_MIN. Rerunning. Not if we set the onlyCancer = FALSE. Also, not if
## we use a sampleEvery = 0.1. Setting it to 0.1 still same issue with
## seed 5613635 or 465554. Again, estimate is less < 1/10000. So a rare
## numerical issue. OK. Decrease sampleEvery further, but increase keepEvery for size.

test_that("mutationPropGrowth no warning with Exp, indiv", {
    ## I once (out of over > 10000) saw it fail. Try to catch it
    
    
    cat("\n a runif is", runif(1), "\n")
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    expect_silent(oiI1 <- oncoSimulIndiv(oi,
                                         model = "Exp",
                                         onlyCancer = FALSE, detectionProb = NA,
                                         sampleEvery = 0.03,
                                         keepEvery = 5,
                                         mutationPropGrowth = TRUE,
                                         seed = NULL))                

})


test_that("mutationPropGrowth no warning with McFl, indiv", {
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
               drvNames = c("d1", "d2", "f1", "f2", "f3"))
              expect_silent(oiI1 <- oncoSimulIndiv(oi,
                                                   model = "McFL",
                                                   mu = 5e-6,
                                                   detectionSize = 1e8, 
                                                   detectionDrivers = 2,
                                                   sampleEvery = 0.025,
                                                   keepEvery = 2,
                                                   onlyCancer = FALSE, detectionProb = NA,
                                                   mutationPropGrowth = TRUE))                
})



cat(paste("\n Ending warning-mutPropGrowth tests", date(), "\n"))



cat(paste("\n Starting warning-mutPropGrowth tests", date(), "\n"))

RNGkind("Mersenne-Twister")

test_that("mutationPropGrowth warning with Bozic, indiv", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = runif(5, 0.01, 0.06),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulIndiv(oi,
                                                    sampleEvery = 0.03,
                                                    keepEvery = 2,
                                                    model = "Bozic",
                                                    mutationPropGrowth = TRUE),
                             "Using fitness Bozic (bozic1) with mutationPropGrowth = TRUE;",
                             fixed = TRUE)                                                     
})


## test_that("mutationPropGrowth warning with Bozic, population", {
##     ## Will not return the warning. Seems a bug in mclapply not respecting
##     ## mc.silent.
##     ## http://stackoverflow.com/questions/21486658/warnings-suppressed-with-mclapply-in-r
##               oi <- allFitnessEffects(orderEffects =
##                c("F > D" = -0.3, "D > F" = 0.4),
##                noIntGenes = rexp(5, 10),
##                           geneToModule =
##                               c("Root" = "Root",
##                                 "F" = "f1, f2, f3",
##                                 "D" = "d1, d2") )
##               expect_warning(oiI1 <- oncoSimulPop(4,
##                                                   oi,
##                                                   model = "Bozic",
##                                                   mutationPropGrowth = TRUE),
##                              "Using fitness Bozic (bozic1) with mutationPropGrowth = TRUE;",
##                              fixed = TRUE)                                                     
##           })

test_that("mutationPropGrowth warning with Bozic, sample", {
    ## Bozic will crash on purpose if any s >= 1. Make sure
    ## that does not happen here
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = runif(5, 0.01, 0.06),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulSample(4,
                                                     oi,
                                                     sampleEvery = 0.03,
                                                     onlyCancer = FALSE,
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
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n the seed is", pseed, "\n")
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2") )
    expect_silent(oiI1 <- oncoSimulIndiv(oi,
                                         model = "Exp",
                                         onlyCancer = TRUE,
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
                                "D" = "d1, d2") )
              expect_silent(oiI1 <- oncoSimulIndiv(oi,
                                                   model = "McFL",
                                                   mu = 5e-6,
                                                   detectionSize = 1e8, 
                                                   detectionDrivers = 2,
                                                   sampleEvery = 0.025,
                                                   keepEvery = 2,
                                                   onlyCancer = TRUE,
                                                   mutationPropGrowth = TRUE))                
})




test_that("mutationPropGrowth OK  with Exp, sample", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_message(oiI1 <- oncoSimulSample(4,
                                                     oi,
                                                     onlyCancer = FALSE,
                                                     model = "Exp",
                                                     sampleEvery = 0.1,
                                                     mutationPropGrowth = TRUE),
                             "Successfully", fixed = TRUE)
})


test_that("mutationPropGrowth OK  with McFL, sample", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_message(oiI1 <- oncoSimulSample(4,
                                                     oi,
                                                     model = "McFL",
                                                     mu = 5e-6,
                                                     detectionSize = 1e8, 
                                                     detectionDrivers = 2,
                                                     sampleEvery = 0.025,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = TRUE),
                             "Successfully", fixed = TRUE)
})

cat(paste("\n Ending warning-mutPropGrowth tests", date(), "\n"))



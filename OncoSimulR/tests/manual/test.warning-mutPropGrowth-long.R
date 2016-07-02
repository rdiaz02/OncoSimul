cat(paste("\n Starting warning-mutPropGrowth-long tests", date(), "\n"))


test_that("mutationPropGrowth OK  with Exp, sample", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2"),
               drvNames = c("d1", "d2", "f1", "f2", "f3"))
              expect_message(oiI1 <- oncoSimulSample(4,
                                                     oi,
                                                     onlyCancer = FALSE, detectionProb = NA,
                                                     model = "Exp",
                                                     sampleEvery = 0.05,
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
                                "D" = "d1, d2"),
              drvNames = c("d1", "d2", "f1", "f2", "f3") )
              expect_message(oiI1 <- oncoSimulSample(4,
                                                     oi,
                                                     model = "McFL",
                                                     mu = 5e-6,
                                                     detectionSize = 1e8, 
                                                     detectionDrivers = 2,
                                                     sampleEvery = 0.025,
                                                     onlyCancer = FALSE, detectionProb = NA,
                                                     mutationPropGrowth = TRUE),
                             "Successfully", fixed = TRUE)
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

cat(paste("\n Ending warning-mutPropGrowth-long tests", date(), "\n"))

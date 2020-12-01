## Test for mutator functionality with FDF
inittime <- Sys.time()
cat(paste("\n Starting test.mutatorFDF.R test at", date()))
date()

## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

date()
test_that("eval fitness and mut OK", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*(f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel")
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  expect_output(ou <- evalGenotypeFitAndMut(genotype = "A",
                                            fitnessEffects = fe,
                                            mutatorEffects = mt,
                                            spPopSizes = c(5000, 2500, 2500, 7500),
                                            verbose = TRUE),
                regexp = "0.5",
                fixed = TRUE)
  
  expect_identical(ou, c(1.5, 1))
  expect_identical(evalGenotypeFitAndMut("B", fe, mt, c(5000,2500,2500,7500)),
                   c(1.5, 10))
  expect_identical(evalGenotypeFitAndMut("A, B", fe, mt, c(5000,2500,2500,7500)),
                   c(5*(5000/(5000+2500+2500+7500)) -
                       0.5*( (2500/(5000+2500+2500+7500)) + (2500/(5000+2500+2500+7500)) ) +
                       15*(7500/(5000+2500+2500+7500)), 10))
  
  
  expect_error(evalGenotypeFitAndMut("A", fe, mt),
               "You have a NULL spPopSizes",
               fixed = TRUE)
  expect_error(evalGenotypeFitAndMut("A", fe, mt, c(5000)),
               "spPopSizes must be as long as number of genotypes",
               fixed = TRUE)

})

test_that("eval mut genotypes", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel")
  
  r2 <- rfitness(2)
  
  fe2 <- allFitnessEffects(genotFitness = r2, 
                           frequencyDependentFitness = FALSE)
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  evge1 <- evalGenotypeFitAndMut(genotype = "Root", fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge2 <- evalGenotypeFitAndMut(genotype = 0, fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge3 <- evalGenotypeFitAndMut(genotype = "A", fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge4 <- evalGenotypeFitAndMut(genotype = 1, fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge5 <- evalGenotypeFitAndMut(genotype = c(1, 2), fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge6 <- evalGenotypeFitAndMut(genotype = "A, B", fitnessEffects = fe,
                                 mutatorEffects = mt,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  expect_equal(evge1, evge2)
  
  expect_equal(evge3, evge4)
  
  expect_equal(evge5, evge6)
  
  expect_error(evalGenotypeFitAndMut(genotype = c(0, 1), fitnessEffects = fe,
                                     mutatorEffects = mt,
                                     spPopSizes = c(5000, 2500, 2500, 7500)), 
               "Genotype cannot contain any 0 if its length > 1")
  
  expect_error(evalGenotypeFitAndMut(genotype = c(1, 3), fitnessEffects = fe,
                                     mutatorEffects = mt,
                                     spPopSizes = c(5000, 2500, 2500, 7500)),
               "Genotype as vector of numbers contains genes not in fitnessEffects/mutatorEffects.")
  
  expect_error(evalGenotypeFitAndMut(genotype = 0, fitnessEffects = fe2,
                                     mutatorEffects = mt), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotypeFitAndMut(genotype = "Root", fitnessEffects = fe2,
                                     mutatorEffects = mt), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotypeFitAndMut(genotype = c(1, 2), fitnessEffects = fe,
                            mutatorEffects = mt), 
               "You have a NULL spPopSizes")
  
  expect_error(evalGenotypeFitAndMut(genotype = "1", fitnessEffects = fe,
                                     mutatorEffects = mt), 
               "You have a NULL spPopSizes")
  
  expect_error(evalGenotypeFitAndMut(genotype = "1", fitnessEffects = fe2), 
               "Genotype contains NAs or genes not in fitnessEffects/mutatorEffects")
})

test_that("eval mut genotypes, echo", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel")
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  expect_output(evalGenotypeFitAndMut("A, B", fe, mt,
                                      spPopSizes = c(5000, 2500, 2500, 7500),
                                      echo = TRUE),
                "Genotype: ", 
                fixed = TRUE)
})

test_that("testing all genes evaluation", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel")
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  genotypes <- c(0, OncoSimulR:::generateAllGenotypes(fitnessEffects = fe, 
                                                      order = FALSE,
                                                      max = 256)$genotNums)
  
  evalGs_one_by_one1 <- sapply(genotypes, function(x) evalGenotypeFitAndMut(x,
                                                                           fitnessEffects = fe,
                                                                           mutatorEffects = mt,
                                                                           spPopSizes = c(5000, 2500, 2500, 7500)))[1, ]
  
  evalGs_all_together1 <- evalAllGenotypesFitAndMut(fitnessEffects = fe,
                                                    mutatorEffects = mt,
                                                    spPopSizes = c(5000, 2500, 2500, 7500),
                                                    )$Fitness
  
  evalGs_one_by_one2 <- sapply(genotypes, function(x) evalGenotypeFitAndMut(x,
                                                                            fitnessEffects = fe,
                                                                            mutatorEffects = mt,
                                                                            spPopSizes = c(5000, 2500, 2500, 7500)))[2, ]
  evalGs_all_together2 <- evalAllGenotypesFitAndMut(fitnessEffects = fe,
                                                    mutatorEffects = mt,
                                                    spPopSizes = c(5000, 2500, 2500, 7500),
                                                    )$MutatorFactor
  
  
  expect_identical(evalGs_one_by_one1, evalGs_all_together1)
  expect_identical(evalGs_one_by_one2, evalGs_all_together2)
})

test_that("expect output oncoSimulIndiv", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel")
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  expect_output(print(oncoSimulIndiv(fe, muEF = mt, sampleEvery = 0.01,
                                     keepEvery = 5)),
                "Individual OncoSimul trajectory",
                fixed = TRUE)
})
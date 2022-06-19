## Test for mutator functionality with FDF
inittime <- Sys.time()
cat(paste("\n Starting test.Z-mutatorFDF at", date(), "\n"))

test_that("Mutator genes missing from fitness", {
  
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  set.seed(1) ## for reproducibility
  
  suppressWarnings(fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel"))
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "C" = 10))
  
  expect_error(oncoSimulIndiv(fe, muEF = mt),
               "Genes in mutatorEffects not present in fitnessEffects",
               fixed = TRUE)
})

test_that("fit. mut. eff. values, long ex",  {
    ## Because of testthat's reluctance to behave sensibly
    silent_expect_true <- function(x) {
        expect_true(suppressWarnings(x))
    }
    
  r1 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(1.5, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  set.seed(1) ## for reproducibility
  
  suppressWarnings(fe <- allFitnessEffects(genotFitness = r1, 
                          frequencyDependentFitness = TRUE, 
                          frequencyType = "rel"))
  
  mt <- allMutatorEffects(epistasis = c("A" = 1, "B" = 10))
  
  silent_expect_true(all.equal(target = evalGenotype("A", fe,
                                              spPopSizes = c(5000, 2500, 2500, 7500)), current = 1.5))
  
  silent_expect_true(all.equal(target = evalGenotype("B", fe,
                                              spPopSizes = c(5000, 2500, 2500, 7500)), current = 1.5))
  
  silent_expect_true(all.equal(target = round(
                                      evalGenotype("A, B", fe,
                                      spPopSizes = c(5000, 2500, 2500, 7500)), 2),
                        current = 7.71))
  
  set.seed(2)
  
  silent_expect_true(all.equal(target = evalGenotype("A", fe,
                                              spPopSizes = c(5000, 2500, 2500, 7500)), current = 1.5))
  
  silent_expect_true(all.equal(target = evalGenotype("B", fe,
                                              spPopSizes = c(5000, 2500, 2500, 7500)), current = 1.5))
  
  silent_expect_true(all.equal(target = round(
                                       evalGenotype("A, B", fe,
                                                    spPopSizes = c(5000, 2500, 2500, 7500)), 2),
                        current = 7.71))
  
  silent_expect_true(all.equal(target = evalGenotypeFitAndMut("A", fe, mt,
                                                       spPopSizes = c(5000, 2500, 2500, 7500)), 
                        current = c(1.5, 1.0)))
  
  silent_expect_true(all.equal(target = evalGenotypeFitAndMut("B", fe, mt,
                                                       spPopSizes = c(5000, 2500, 2500, 7500)), 
                        current = c(1.5, 10)))
  
  silent_expect_true(all.equal(target = round(
                                       evalGenotypeFitAndMut("A, B", fe, mt,
                                                             spPopSizes = c(5000, 2500, 2500, 7500)), 2),
                        current = c(7.71, 10)))
})

set.seed(NULL)

cat(paste("\n Ending test.Z-mutatorFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

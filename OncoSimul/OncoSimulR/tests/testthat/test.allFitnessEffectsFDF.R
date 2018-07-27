inittime <- Sys.time()
cat(paste("\n Starting test.allFitnessEffectsFDF at", date(), "\n"))

test_that("Errors expectations", {
  r1 <- rfitness(2)
  r2 <- data.frame(r1)
  r3 <- r2
  r3[, "Fitness"] <- c("1", "2", "3", "4")
  r4 <- r2
  r4[, "Fitness"] <- c("F_", "F_1", "F_2", "F_1_2")
  r5 <- r2
  r5[, "Fitness"] <- c("f_", "f_1", "f_2", "f_1_2")
  r6 <- NULL
  r7 <- data.frame()
  
  expect_error(allFitnessEffects(genotFitness = r1, 
                                   frequencyDependentFitness = TRUE), 
               "Input must inherit from data.frame.")
  
  expect_error(allFitnessEffects(genotFitness = r2, 
                                 frequencyDependentFitness = TRUE), 
               "All elements in last column must be character.")
  
  expect_error(allFitnessEffects(genotFitness = r3, 
                                 frequencyDependentFitness = TRUE), 
               "There are some errors in fitness column")
  
  expect_error(allFitnessEffects(genotFitness = r4, 
                                 frequencyDependentFitness = TRUE), 
               "There are some errors in fitness column")
  
  expect_error(allFitnessEffects(genotFitness = r5, 
                                 frequencyDependentFitness = FALSE), 
               "A genotype fitness matrix/data.frame must be numeric.")
  
  
  expect_error(allFitnessEffects(genotFitness = r6, 
                                 frequencyDependentFitness = TRUE), 
               "You have a null genotFitness in a frequency dependent fitness situation.")
  
  expect_error(allFitnessEffects(genotFitness = r7, 
                                 frequencyDependentFitness = TRUE), 
               "You have an empty data.frame")
  
  
})

test_that("testing output", {
  
  r1 <- data.frame(rfitness(3))
  r1[, "Fitness"] <- c("max(1, f_)",
                      "max(1, f_1)", 
                      "max(1, f_2)", 
                      "max(1, f_3)",
                      "100*f_1_2 + max(max(1, f_1), max(1, f_2))",
                      "100*f_1_3 + max(max(1, f_1), max(1, f_3))",
                      "100*f_2_3 + max(max(1, f_2), max(1, f_3))",
                      "200*f_1_2_3 + 500*(f_1_2 + f_1_3 + f_2_3)") 
  
  afe <- allFitnessEffects(genotFitness = r1, 
                           frequencyDependentFitness = TRUE)
  
  expect_true(afe$frequencyDependentFitness)
  
  lapply(afe[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe$gMOneToOne)
  
  lapply(afe[c(10:13, 19)],  function(x){
    expect_null(x)
  })
  
  expect_equivalent(afe$geneToModule, "Root")
  
  expect_identical(afe$fitnessLandscapeVariables, 
                   OncoSimulR:::fVariables(ncol(afe$fitnessLandscape) - 1))
  
  lapply(afe[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  expect_length(afe$drv, 0)
  
  expect_equal(class(afe$drv), "integer")
  
})

cat(paste("\n Ending test.allFitnessEffectsFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

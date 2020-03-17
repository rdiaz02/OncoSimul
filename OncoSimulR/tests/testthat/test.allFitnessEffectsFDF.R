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
  r8 <- r2
  r8[, "Fitness"] <- c("n_", "n_1", "n_2", "n_1_2")
  
  expect_error(allFitnessEffects(genotFitness = r1, 
                                 frequencyDependentFitness = TRUE, 
                                 frequencyType = "rel"), 
               "Input must inherit from data.frame.")
  
  expect_error(allFitnessEffects(genotFitness = r2, 
                                 frequencyDependentFitness = TRUE, 
                                 frequencyType = "rel"), 
               "All elements in last column must be character.")
  
  expect_error(allFitnessEffects(genotFitness = r3, 
                                 frequencyDependentFitness = TRUE, 
                                 frequencyType = "rel"), 
               "There are some errors in fitness column")
  
  expect_error(allFitnessEffects(genotFitness = r4, 
                                 frequencyDependentFitness = TRUE, 
                                 frequencyType = "rel"), 
               "There are some errors in fitness column")
  
  expect_error(allFitnessEffects(genotFitness = r5, 
                                 frequencyDependentFitness = FALSE), 
               "A genotype fitness matrix/data.frame must be numeric.")
  
  
  expect_error(allFitnessEffects(genotFitness = r6, 
                                 frequencyDependentFitness = TRUE,
                                 frequencyType = "rel"),  
               "You have a null genotFitness in a frequency dependent fitness situation.")
  
  expect_error(allFitnessEffects(genotFitness = r7, 
                                 frequencyDependentFitness = TRUE,
                                 frequencyType = "rel"),  
               "You have an empty data.frame")
  
  expect_error(allFitnessEffects(genotFitness = r8, 
                                 frequencyDependentFitness = TRUE,
                                 frequencyType = "rel"),  
               "There are some errors in fitness column")
  
  expect_error(allFitnessEffects(genotFitness = r8, 
                                 frequencyDependentFitness = TRUE),  
               "frequencyType must be 'abs' \\(absolute\\) or 'rel' \\(relative\\).")
  
  
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
  
  r2 <- rfitness(2)
  
  r3 <- data.frame(A = c(0, 1, 0, 1),
                   B = c(0, 0, 1, 1),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"))
  
  r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"))
  
  r5 <- data.frame(A = c(0, 1, 0, 1),
                   B = c(0, 0, 1, 1),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_A))",
                               "max(2, 3*(f_ + f_B))",
                               "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"))
  
  r6 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_A))",
                               "max(2, 3*(f_ + f_B))",
                               "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"))
  
  r7 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*n_)",
                               "max(1.5, 3*(n_ + n_A))",
                               "max(2, 3*(n_ + n_B))",
                               "max(2, 5*n_ - 0.5*( n_A + n_B) + 15*n_A_B)"))

  afe1 <- allFitnessEffects(genotFitness = r1, 
                            frequencyDependentFitness = TRUE,
                            frequencyType = "rel")
  
  afe2 <- allFitnessEffects(genotFitness = r7, 
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")
  
  
  expect_true(afe1$frequencyDependentFitness)
  
  expect_true(afe2$frequencyDependentFitness)
  
  lapply(afe1[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  lapply(afe2[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe1$gMOneToOne)
  
  expect_true(afe2$gMOneToOne)
  
  lapply(afe1[c(10:13, 20)],  function(x){
    expect_null(x)
  })
  
  lapply(afe2[c(10:13, 20)],  function(x){
    expect_null(x)
  })
  
  expect_equivalent(afe1$geneToModule, "Root")
  
  expect_equivalent(afe2$geneToModule, "Root")
  
  expect_identical(afe1$fitnessLandscapeVariables, 
                   OncoSimulR:::fVariablesN(ncol(afe1$fitnessLandscape) - 1, 
                                            frequencyType = "rel"))
  
  expect_identical(afe2$fitnessLandscapeVariables, 
                   OncoSimulR:::fVariablesN(ncol(afe2$fitnessLandscape) - 1, 
                                            frequencyType = "abs"))
  
  lapply(afe1[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  lapply(afe2[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  expect_length(afe1$drv, 0)
  
  expect_length(afe2$drv, 0)
  
  expect_equal(class(afe1$drv), "integer")
  
  expect_equal(class(afe2$drv), "integer")
  
  expect_warning(allFitnessEffects(genotFitness = r2, 
                                   frequencyDependentFitness = FALSE, 
                                   frequencyType = "rel",
                                   spPopSizes = c(500, 
                                                  250, 
                                                  250, 
                                                  250, 
                                                  300,
                                                  300,
                                                  300,
                                                  450)) , 
                 "spPopSizes will be considered NULL if frequencyDependentFitness = FALSE")
  if(as.character(version$major) < 4) {
  expect_warning(allFitnessEffects(genotFitness = r3, 
                                   frequencyDependentFitness = TRUE, 
                                   frequencyType = "rel"), 
                 "Last column of genotype fitness is a factor. Converting to character.")
  }

  expect_identical(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
               allFitnessEffects(genotFitness = r4, 
                                 frequencyDependentFitness = TRUE, 
                                 frequencyType = "rel"))
  
  expect_identical(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
                   allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"))
  
  expect_identical(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
                   allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"))
  
  expect_identical(allFitnessEffects(genotFitness = r4, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
                   allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"))
  
  expect_identical(allFitnessEffects(genotFitness = r4, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
                   allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"))
  
  expect_identical(allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"), 
                   allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE, 
                                     frequencyType = "rel"))
})

cat(paste("\n Ending test.allFitnessEffectsFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

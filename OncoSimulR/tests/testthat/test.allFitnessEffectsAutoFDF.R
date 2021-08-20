inittime <- Sys.time()
cat(paste("\n Starting test.allFitnessEffectsFDF at", date(), "\n"))

test_that("testing output when frequencyType = 'auto'", {
  
  r1 <- data.frame(rfitness(3))
  
  colnames(r1)[which(colnames(r1) == "Birth")] <- "Fitness"
  
  r1[, "Fitness"] <- c("max(1, f_)",
                       "max(1, f_1)", 
                       "max(1, f_2)", 
                       "max(1, f_3)",
                       "100*f_1_2 + max(max(1, f_1), max(1, f_2))",
                       "100*f_1_3 + max(max(1, f_1), max(1, f_3))",
                       "100*f_2_3 + max(max(1, f_2), max(1, f_3))",
                       "200*f_1_2_3 + 500*(f_1_2 + f_1_3 + f_2_3)") 
  
  r2 <- rfitness(2)
  colnames(r2)[which(colnames(r2) == "Birth")] <- "Fitness"
  
  r3 <- data.frame(A = c(0, 1, 0, 1),
                   B = c(0, 0, 1, 1),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
                   stringsAsFactors = FALSE)
  
  r5 <- data.frame(A = c(0, 1, 0, 1),
                   B = c(0, 0, 1, 1),
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_A))",
                               "max(2, 3*(f_ + f_B))",
                               "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                   stringsAsFactors = FALSE)
  
  r6 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_A))",
                               "max(2, 3*(f_ + f_B))",
                               "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                   stringsAsFactors = FALSE)
  
  r7 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("max(3, 2*n_)",
                               "max(1.5, 3*(n_ + n_A))",
                               "max(2, 3*(n_ + n_B))",
                               "max(2, 5*n_ - 0.5*( n_A + n_B) + 15*n_A_B)"),
                   stringsAsFactors = FALSE)
  
  
  suppressWarnings(afe1 <- allFitnessEffects(genotFitness = r1, 
                            frequencyDependentFitness = TRUE))
  
  suppressWarnings(afe2 <- allFitnessEffects(genotFitness = r7, 
                            frequencyDependentFitness = TRUE))
  
  expect_message(suppressWarnings(allFitnessEffects(genotFitness = r1,
                                   frequencyDependentFitness = TRUE)),
                 "frequencyType set to 'auto'")
  
  expect_message(suppressWarnings(allFitnessEffects(genotFitness = r7, 
                                   frequencyDependentFitness = TRUE)),
                 "frequencyType set to 'auto'")
  
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE)),
                   suppressWarnings(allFitnessEffects(genotFitness = r4,
                                     frequencyDependentFitness = TRUE)))

  ## 20 is where we are storing full_FDF_spec for now. Very brittle
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]),
                   suppressWarnings(allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]))
  
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r3, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]),
                   suppressWarnings(allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]))
  
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r4, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]),
                   suppressWarnings(allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]))
  
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r4, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]),
                   suppressWarnings(allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE)[-c(14, 19)]))
  
  expect_identical(suppressWarnings(allFitnessEffects(genotFitness = r5, 
                                     frequencyDependentFitness = TRUE)),
                   suppressWarnings(allFitnessEffects(genotFitness = r6, 
                                     frequencyDependentFitness = TRUE)))

})

cat(paste("\n Ending test.allFitnessEffectsFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

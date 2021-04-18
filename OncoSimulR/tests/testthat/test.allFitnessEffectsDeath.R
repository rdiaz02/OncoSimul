inittime <- Sys.time()
cat(paste("\n Starting test.allFitnessEffectsDeath at", date(), "\n"))

test_that("Errors expectations", {
  r1 <- rfitness(2)
  r2 <- data.frame(r1)
  r2[, "Death"] <- c(1, 2, 3, 4)
  # No FDB, no FDD
  r3 <- r2
  r3[, "Birth"] <- c(1, 2, 3, 4)
  r3[, "Death"] <- c(1, 2, 3, 4)
  # Expect error
  r4 <- r2
  r4[, "Birth"] <- c(1, 2, 3, 4)
  r4[, "Death"] <- c("F_", "F_1", "F_2", "F_1_2")
  # Expect error
  r5 <- r2
  r5[, "Birth"] <- c(1, 2, 3, 4)
  r5[, "Death"] <- c("N_", "N_1", "N_2", "N_1_2")
  # FDB, no FDD, rel freq
  r6 <- r2
  r6[, "Birth"] <- c("f_", "f_1", "f_2", "f_1_2")
  r6[, "Death"] <- c(1, 2, 3, 4)
  # No FDB, FDD, rel freq
  r7 <- r2
  r7[, "Birth"] <- c(1, 2, 3, 4)
  r7[, "Death"] <- c("f_", "f_1", "f_2", "f_1_2")

  # FDB, FDD, abs freq
  r8 <- r2
  r8[, "Birth"] <- c("n_", "n_1", "n_2", "n_1_2")
  r8[, "Death"] <- c("n_", "n_1", "n_2", "n_1_2")
  
  # Expect error
  r9 <- r2
  r9[, "Birth"] <- c("f_", "f_1", "f_2", "f_1_2")
  r9[, "Death"] <- c("n_", "n_1", "n_2", "n_1_2")
  
  expect_error(allFitnessEffects(genotFitness = r1, 
                                 frequencyDependentBirth = TRUE, 
                                 frequencyType = "rel"), 
               "Input must inherit from data.frame.")
  
  expect_error(allFitnessEffects(genotFitness = r1, 
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"), 
               "genotFitness: if genotype is specified, it must be data frame")
  
  expect_error(allFitnessEffects(genotFitness = r2, 
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"), 
               "All elements in last column must be character.")
  
  expect_error(allFitnessEffects(genotFitness = r2, 
								 frequencyDependentBirth = TRUE,
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"), 
               "All elements in last column must be character.")
  
  expect_error(allFitnessEffects(genotFitness = r3, 
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,								 
                                 frequencyType = "rel"), 
                "All elements in last column must be character.")
  
  # Note that if Birth has 0 or 1 as value in each row, this error
  # will not trigger and the function will work. There will be
  # 3 genes (A, B and Birth) and the Birth Rate will be the one
  # in Death column. Maybe this is confusing and we need to
  # warn the user when the naming of the last column is Death, but
  # deathSpec = FALSE.
  expect_error(allFitnessEffects(genotFitness = r3), 
               "First ncol - 1 entries not in \\{0, 1\\}.")
			   
  expect_error(allFitnessEffects(genotFitness = r4, 
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"), 
               "There are some errors in death column")
			   
	expect_error(allFitnessEffects(genotFitness = r5, 
					 frequencyDependentDeath = TRUE,
					 deathSpec = TRUE,
					 frequencyType = "rel"), 
	"There are some errors in death column")
  
  expect_error(allFitnessEffects(genotFitness = r6, 
                                 frequencyDependentBirth = FALSE,
								 deathSpec = TRUE), 
               "A genotype fitness matrix/data.frame must be numeric.")
  
  expect_error(allFitnessEffects(genotFitness = r7, 
                                 frequencyDependentDeath = FALSE,
								 deathSpec=TRUE), 
               "A genotype fitness matrix/data.frame must be numeric.")
  
  expect_error(allFitnessEffects(genotFitness = r7, 
                                 frequencyDependentBirth = TRUE,
								 deathSpec = TRUE), 
               "All columns except from the second last must be numeric.")
  
  expect_error(allFitnessEffects(genotFitness = r8,
								 frequencyDependentBirth = TRUE,
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"),  
               "There are some errors in birth column")
			   
  expect_error(allFitnessEffects(genotFitness = r9,
								 frequencyDependentBirth = TRUE,
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE,
                                 frequencyType = "rel"),  
               "There are some errors in death column")
  
  expect_error(allFitnessEffects(genotFitness = r9,
								 frequencyDependentBirth = TRUE,
                                 frequencyDependentDeath = TRUE,
								 deathSpec = TRUE),
               "There are some errors in death column")

  
  
})

cat(paste("\n Ending test.allFitnessEffectsDeath at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

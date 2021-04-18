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

test_that("testing output", {
  
  r1 <- data.frame(rfitness(3))
  
  r1[, "Birth"] <- c("max(1, f_)",
                       "max(1, f_1)", 
                       "max(1, f_2)", 
                       "max(1, f_3)",
                       "100*f_1_2 + max(max(1, f_1), max(1, f_2))",
                       "100*f_1_3 + max(max(1, f_1), max(1, f_3))",
                       "100*f_2_3 + max(max(1, f_2), max(1, f_3))",
                       "200*f_1_2_3 + 500*(f_1_2 + f_1_3 + f_2_3)") 
  
  r2 <- rfitness(2)
  
  r3 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Birth = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
				   Death = c(1, 2, 3, 4))
  
  r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Birth = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"),
				   Death = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_1))",
                               "max(2, 3*(f_ + f_2))",
                               "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)"))
  
  r5 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
				   Birth = c(1, 2, 3, 4),
                   Death = c(1, 2, 3, 4))
							   
  r6 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
				   Birth = c(1, 2, 3, 4),
                   Death = c("max(3, 2*f_)",
                               "max(1.5, 3*(f_ + f_A))",
                               "max(2, 3*(f_ + f_B))",
                               "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"))
  

  afe1 <- allFitnessEffects(genotFitness = r3, 
                            frequencyDependentBirth = TRUE,
							deathSpec = TRUE,
                            frequencyType = "rel")
  
  
  expect_true(afe1$frequencyDependentBirth)
  expect_false(afe1$frequencyDependentDeath)
  expect_true(afe1$deathSpec)
  
  lapply(afe1[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe1$gMOneToOne)
  
  lapply(afe1[c(10:13)],  function(x){
    expect_null(x)
  })

  expect_equivalent(afe1$geneToModule, "Root")
  
  expect_equivalent(afe1$full_FDF_spec[, "Genotype_as_fvarsb"],
                   OncoSimulR:::fVariablesN(ncol(afe1$fitnessLandscape) - 2, 
                                            frequencyType = "rel"))
  
  lapply(afe1[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  
  expect_length(afe1$drv, 0)
  
  expect_equal(class(afe1$drv), "integer")
  
  afe2 <- allFitnessEffects(genotFitness = r4, 
							frequencyDependentBirth = TRUE,
							frequencyDependentDeath = TRUE,
							deathSpec = TRUE,
							frequencyType = "rel")
  
  
  expect_true(afe2$frequencyDependentBirth)
  expect_true(afe2$frequencyDependentDeath)
  expect_true(afe2$deathSpec)
  
  lapply(afe2[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe2$gMOneToOne)
  
  lapply(afe2[c(10:13)],  function(x){
    expect_null(x)
  })

  expect_equivalent(afe2$geneToModule, "Root")
  
  expect_equivalent(afe2$full_FDF_spec[, "Genotype_as_fvarsb"],
                   OncoSimulR:::fVariablesN(ncol(afe1$fitnessLandscape) - 2, 
                                            frequencyType = "rel"))
											
  expect_equivalent(afe2$full_FDF_spec[, "Genotype_as_fvarsd"],
                   OncoSimulR:::fVariablesN(ncol(afe1$fitnessLandscape) - 2, 
                                            frequencyType = "rel"))
  
  lapply(afe2[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  expect_length(afe2$drv, 0)
  
  expect_equal(class(afe2$drv), "integer")
 
  expect_warning(afe3 <- allFitnessEffects(genotFitness = r5, 
							frequencyDependentBirth = FALSE,
							frequencyDependentDeath = FALSE,
							deathSpec = TRUE,
							frequencyType = "rel"), "frequencyType set to NA")
	
  
  expect_false(afe3$frequencyDependentBirth)
  expect_false(afe3$frequencyDependentDeath)
  expect_true(afe3$deathSpec)
  
  lapply(afe3[c(1:3)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe3$gMOneToOne)
  
  lapply(afe3[c(8, 10:13)],  function(x){
    expect_null(x)
  })

  expect_equivalent(afe3$geneToModule, "Root")
  
  lapply(afe3[c(4, 5, 15:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  expect_length(afe3$fitnessLandscapeVariables, 0)
  expect_equal(class(afe3$fitnessLandscapeVariables), "character")
  
  expect_length(afe3$drv, 0)
  expect_equal(class(afe3$drv), "integer")
  
  afe4 <- allFitnessEffects(genotFitness = r6, 
                            frequencyDependentDeath = TRUE,
							deathSpec = TRUE,
                            frequencyType = "rel")
  
  
  expect_false(afe4$frequencyDependentBirth)
  expect_true(afe4$frequencyDependentDeath)
  expect_true(afe4$deathSpec)
  
  lapply(afe4[c(1:3, 8)], function(x){
    expect_equal(length(x), 0)
    
    expect_equal(class(x), "list")
  })
  
  expect_true(afe4$gMOneToOne)
  
  lapply(afe4[c(10:13)],  function(x){
    expect_null(x)
  })

  expect_equivalent(afe4$geneToModule, "Root")
  
  expect_equivalent(afe4$full_FDF_spec[, "Genotype_as_fvarsd"],
                   OncoSimulR:::fVariablesN(ncol(afe4$fitnessLandscape) - 2, 
                                            frequencyType = "rel"))
  
  lapply(afe4[c(4, 5, 14:16)], function(x){
    expect_equal(class(x), "data.frame")
  })
  
  
  expect_length(afe4$drv, 0)
  
  expect_equal(class(afe4$drv), "integer")

})

cat(paste("\n Ending test.allFitnessEffectsDeath at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

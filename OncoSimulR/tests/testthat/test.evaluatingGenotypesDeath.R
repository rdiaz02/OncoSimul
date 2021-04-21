inittime <- Sys.time()
cat(paste("\n Starting test.evaluatingGenotypesDeath at", date(), "\n"))

test_that("testing single gene evaluation", {
  
  r1 <- data.frame(rfitness(2))
  
  r1[, "Birth"] <- c("max(3, 2*f_)",
                      "max(1.5, 3*(f_ + f_1))",
                      "max(1.5, 3*(f_ + f_2))",
                      "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)")
  
  r1[, "Death"] <- c(2, 3, 3, 4)
  
  r2 <- rfitness(2)
  
  r2 <- cbind(r2, Death=c(2, 3, 3, 4))
  
  r3 <- data.frame(rfitness(2))
			
  r3[, "Birth"] <- c(2, 3, 3, 4)
  
  r3[, "Death"] <- c("max(3, 2*f_)",
                      "max(1.5, 3*(f_ + f_1))",
                      "max(1.5, 3*(f_ + f_2))",
                      "max(2, 5*f_ - 0.5*( f_1 + f_2) + 15*f_1_2)")
  
  afe1 <- allFitnessEffects(genotFitness = r1, 
                            frequencyDependentBirth = TRUE,
							deathSpec = TRUE,
			                frequencyType = "rel")
  
  suppressWarnings(afe2 <- allFitnessEffects(genotFitness = r2, 
                            frequencyDependentBirth = FALSE,
							deathSpec=TRUE))
  
  afe3 <- allFitnessEffects(genotFitness = r3, 
                            frequencyDependentBirth = FALSE,
							frequencyDependentDeath = TRUE,
							deathSpec=TRUE,
			    frequencyType = "rel")
				
  suppressWarnings({
  evge1 <- evalGenotype(genotype = "Root", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge2 <- evalGenotype(genotype = 0, fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge3 <- evalGenotype(genotype = "A", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge4 <- evalGenotype(genotype = 1, fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge5 <- evalGenotype(genotype = c(1, 2), fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge6 <- evalGenotype(genotype = "A, B", fitnessEffects = afe1,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  })
  
  expect_equal(evge1, evge2)
  
  expect_equal(evge3, evge4)
  
  expect_equal(evge5, evge6)
  
  suppressWarnings({
  expect_error(evalGenotype(genotype = c(0, 1), fitnessEffects = afe1,
                            spPopSizes = c(5000, 2500, 2500, 7500)), 
               "Genotype cannot contain any 0 if its length > 1")
  
  expect_error(evalGenotype(genotype = c(1, 3), fitnessEffects = afe1,
                            spPopSizes = c(5000, 2500, 2500, 7500)),
               "Genotype as vector of numbers contains genes not in fitnessEffects/mutatorEffects.")
			   
  expect_error(evalGenotype(genotype = 0, fitnessEffects = afe2), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotype(genotype = "Root", fitnessEffects = afe2), 
               "Genotype cannot be 0.")
  
  expect_error(evalGenotype(genotype = "1", fitnessEffects = afe1), 
               "You have a NULL spPopSizes")
  
  expect_error(evalGenotype(genotype = "1", fitnessEffects = afe2), 
               "Genotype contains NAs or genes not in fitnessEffects/mutatorEffects")})
  
   suppressWarnings({
  evge1 <- evalGenotype(genotype = "Root", fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge2 <- evalGenotype(genotype = 0, fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge3 <- evalGenotype(genotype = "A", fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge4 <- evalGenotype(genotype = 1, fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge5 <- evalGenotype(genotype = c(1, 2), fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  
  evge6 <- evalGenotype(genotype = "A, B", fitnessEffects = afe3,
                        spPopSizes = c(5000, 2500, 2500, 7500))
  })
  
  expect_equal(evge1, evge2)
  
  expect_equal(evge3, evge4)
  
  expect_equal(evge5, evge6)
})

test_that("testing all genes evaluation", {
  
  r <- data.frame(rfitness(3))
  
  r1 <- r
  
  r1[, "Birth"] <- c("max(2, log(1 + f_))",
                      "3",
                      "3",
                      "3", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
                      "max(5, (1 + f_1_2 + f_2_3 + f_1_3)^2)")
					  
  r1[, "Death"] <- c("max(2, log(1 + f_))",
					  "3",
					  "3",
					  "3", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(5, (1 + f_1_2 + f_2_3 + f_1_3)^2)")
					  
  r2 <- r
  
  r2[, "Birth"] <- c(1, 2, 3, 4, 5, 6, 7, 8)
					  
  r2[, "Death"] <- c("max(2, log(1 + f_))",
					  "3",
					  "3",
					  "3", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(5, (1 + f_1_2 + f_2_3 + f_1_3)^2)")
	
  r3 <- r
  
  r3[, "Birth"] <- c("max(2, log(1 + f_))",
					  "3",
					  "3",
					  "3", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(3, log(1 + f_ + f_1 + f_2 + f_3))", 
					  "max(5, (1 + f_1_2 + f_2_3 + f_1_3)^2)")
					  
  r3[, "Death"] <- c(1, 2, 3, 4, 5, 6, 7, 8)
  
  afe1 <- allFitnessEffects(genotFitness = r1, 
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
			                     frequencyType = "rel")
                           #spPopSizes = c(500, 
                            #              250, 
                             #             250, 
                              #            250, 
                               #           300,
                                #          300,
                                 #         300,
                                  #        450))
  
  genotypes <- c(0, OncoSimulR:::generateAllGenotypes(fitnessEffects = afe1, 
                                                      order = FALSE, 
                                                      max = 256)$genotNums)

  suppressWarnings({
  evalGs_one_by_one <- sapply(genotypes, function(x) evalGenotype(x, afe1,
                                                                  spPopSizes = c(500,
                                                                                 250, 
                                                                                 250, 
                                                                                 250, 
                                                                                 300,
                                                                                 300,
                                                                                 300,
                                                                                 450)))
  
  aux <- evalAllGenotypes(afe1, spPopSizes = c(500,
												 250,
												 250,
												 250,
												 300,
												 300,
												 300,
												 450))
  evalGs_all_together <- rbind(aux$Birth, aux$Death)

  })
  
  expect_identical(evalGs_one_by_one, evalGs_all_together)
  
  afe2 <- allFitnessEffects(genotFitness = r2, 
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
			                     frequencyType = "rel")
                           #spPopSizes = c(500, 
                            #              250, 
                             #             250, 
                              #            250, 
                               #           300,
                                #          300,
                                 #         300,
                                  #        450))
  
  genotypes <- c(0, OncoSimulR:::generateAllGenotypes(fitnessEffects = afe2, 
                                                      order = FALSE, 
                                                      max = 256)$genotNums)

  suppressWarnings({
  evalGs_one_by_one <- sapply(genotypes, function(x) evalGenotype(x, afe2,
                                                                  spPopSizes = c(500,
                                                                                 250, 
                                                                                 250, 
                                                                                 250, 
                                                                                 300,
                                                                                 300,
                                                                                 300,
                                                                                 450)))
  
  aux <- evalAllGenotypes(afe2, spPopSizes = c(500,
												 250,
												 250,
												 250,
												 300,
												 300,
												 300,
												 450))
  evalGs_all_together <- rbind(aux$Birth, aux$Death)

  })
  
  expect_identical(evalGs_one_by_one, evalGs_all_together)
  
  afe3 <- allFitnessEffects(genotFitness = r3, 
                           frequencyDependentBirth = TRUE,
						   deathSpec = TRUE,
			                     frequencyType = "rel")
                           #spPopSizes = c(500, 
                            #              250, 
                             #             250, 
                              #            250, 
                               #           300,
                                #          300,
                                 #         300,
                                  #        450))
  
  genotypes <- c(0, OncoSimulR:::generateAllGenotypes(fitnessEffects = afe3, 
                                                      order = FALSE, 
                                                      max = 256)$genotNums)

  suppressWarnings({
  evalGs_one_by_one <- sapply(genotypes, function(x) evalGenotype(x, afe3,
                                                                  spPopSizes = c(500,
                                                                                 250, 
                                                                                 250, 
                                                                                 250, 
                                                                                 300,
                                                                                 300,
                                                                                 300,
                                                                                 450)))
  
  aux <- evalAllGenotypes(afe3, spPopSizes = c(500,
												 250,
												 250,
												 250,
												 300,
												 300,
												 300,
												 450))
  evalGs_all_together <- rbind(aux$Birth, aux$Death)

  })
  
  expect_identical(evalGs_one_by_one, evalGs_all_together)
  
})


cat(paste("\n Ending test.evaluatingGenotypesDeath at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

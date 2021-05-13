inittime <- Sys.time()
cat(paste("\n Starting test.Z-oncoSimulIndivConstant at", date(), "\n"))

test_that("testing models", {
	
	r <- data.frame(rfitness(2))
	
	r[, "Birth"] <- c("f_", "f_1", "f_2", "f_1_2")
	
	afe1 <- allFitnessEffects(genotFitness = r, 
                            frequencyDependentBirth = TRUE,
						    deathSpec = FALSE)
							
	set.seed(1)
	
	null <- capture.output(osi1 <- oncoSimulIndiv(afe1, 
						model = "Const", 
						onlyCancer = FALSE, 
						finalTime = 100, 
						verbosity = 0, 
						mu = 1e-6,
						initSize = c(1000, 1000),
						initMutant = c("WT", "A"),
						keepPhylog = FALSE,
						seed = NULL, 
						errorHitMaxTries = TRUE, 
						errorHitWallTime = TRUE))
	
	expect_true(osi1$TotalPopSize >= 1800)
	expect_true(osi1$TotalPopSize <= 2200)	
	
	r[, "Death"] <- c("f_1", "f_1 + f_", "f_2+f_1_2", "f_2")
	
	afe2 <- allFitnessEffects(genotFitness = r, 
                            frequencyDependentBirth = TRUE,
							frequencyDependentDeath = TRUE,
						    deathSpec = TRUE)
							
	null <- capture.output(osi2<- oncoSimulIndiv(afe2, 
						model = "Const", 
						onlyCancer = FALSE, 
						finalTime = 100, 
						verbosity = 0, 
						mu = 1e-6,
						initSize = c(1000, 1000),
						initMutant = c("WT", "A"),
						keepPhylog = FALSE,
						seed = NULL, 
						errorHitMaxTries = TRUE, 
						errorHitWallTime = TRUE))
	
	expect_true(osi2$TotalPopSize >= 1800)
	expect_true(osi2$TotalPopSize <= 2200)	
})

set.seed(NULL)

cat(paste("\n Ending test.Z-oncoSimulIndivConstant at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
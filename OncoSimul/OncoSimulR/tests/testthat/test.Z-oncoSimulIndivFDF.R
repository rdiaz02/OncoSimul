inittime <- Sys.time()
cat(paste("\n Starting test.Z-oncoSimulIndivFDF at", date(), "\n"))


test_that("testing output classes", {
  
  r <- data.frame(rfitness(2))
  
  r[, "Fitness"] <- c("f_ - f_1 - f_2 - f_1_2", 
                      "max(100*f_1, 10)", 
                      "max(100*f_2, 10)", 
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")
  
  
  afe <- allFitnessEffects(genotFitness = r, 
                           frequencyDependentFitness = TRUE, 
                           frequencyType = "rel")
  
  set.seed(1)
  
  osi <- oncoSimulIndiv(afe, 
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 5000, 
                        verbosity = 0, 
                        mu = 1e-6,
                        initSize = 500, 
                        keepPhylog = FALSE,
                        seed = NULL, 
                        errorHitMaxTries = TRUE, 
                        errorHitWallTime = TRUE)
  
  expect_identical(class(r), "data.frame")
  
  expect_identical(class(afe), "fitnessEffects")
  
  
  expect_identical(class(osi), c("oncosimul", "oncosimul2"))
  
  expect_identical(class(osi$Genotypes), "matrix")
  
})

test_that("testing performance", {
  
  r <- data.frame(rfitness(2))
  
  r[, "Fitness"] <- c("10*f_", 
                      "10*f_1", 
                      "50*f_2", 
                      "200*(f_1 + f_2) + 50*f_1_2")
  
  
  afe <- allFitnessEffects(genotFitness = r, 
                           frequencyDependentFitness = TRUE, 
                           frequencyType = "rel")
  
  set.seed(1)
  
  osi <- oncoSimulIndiv(afe, 
                        model = "Exp", 
                        onlyCancer = FALSE, 
                        finalTime = 5000, 
                        verbosity = 0, 
                        mu = 1e-6,
                        initSize = 500, 
                        keepPhylog = FALSE,
                        seed = NULL, 
                        errorHitMaxTries = TRUE, 
                        errorHitWallTime = TRUE)
  
  expect_output(print(osi),
                "Individual OncoSimul trajectory", fixed = TRUE)
  
  expect_identical(osi$NumClones, 3)
  
  expect_true(osi$TotalPopSize >2e8)
  
  expect_true(osi$geneNames[2] %in% LETTERS)
  
  expect_false(osi$FinalTime < 1.4)
  
  expect_equal(osi$NumIter, 458)
  
})

test_that("testing Bozic failure", {
  r <- data.frame(rfitness(2))
  
  r[, "Fitness"] <- c("10*f_", 
                      "10*f_1", 
                      "50*f_2", 
                      "200*(f_1 + f_2) + 50*f_1_2")
  
  
  afe <- allFitnessEffects(genotFitness = r, 
                           frequencyDependentFitness = TRUE, 
                           frequencyType = "rel")
  
  set.seed(1)
  
  st <- capture.output(osi <- oncoSimulIndiv(afe, 
                                             model = "Bozic", 
                                             onlyCancer = FALSE, 
                                             finalTime = 5000, 
                                             verbosity = 0, 
                                             mu = 1e-6,
                                             initSize = 500, 
                                             keepPhylog = FALSE,
                                             seed = NULL, 
                                             errorHitMaxTries = TRUE, 
                                             errorHitWallTime = TRUE))
  expect_true(st[22] == " Unrecoverable exception: Algo 2: retval not finite. Aborting. ")

})


set.seed(NULL)

cat(paste("\n Ending test.Z-oncoSimulIndivFDF at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

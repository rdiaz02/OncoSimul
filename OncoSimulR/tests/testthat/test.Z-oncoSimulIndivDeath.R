inittime <- Sys.time()
cat(paste("\n Starting test.Z-oncoSimulIndivDeath at", date(), "\n"))


test_that("testing model and afe compatibility", {

  r <- data.frame(rfitness(2))

  r[, "Birth"] <- c("f_ - f_1 - f_2 - f_1_2",
                      "max(100*f_1, 10)",
                      "max(100*f_2, 10)",
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")


  afe <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
                           frequencyType = "rel")

  expect_error(oncoSimulIndiv(afe,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 20,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 5000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE),
	"To use Arb model specify both birth and death in fitness effects.")

  r[, "Death"] <- c("f_ - f_1 - f_2 - f_1_2",
                      "max(100*f_1, 10)",
                      "max(100*f_2, 10)",
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")

  afe <- allFitnessEffects(genotFitness = r,
					   frequencyDependentBirth = TRUE,
					   frequencyDependentDeath = TRUE,
					   deathSpec = TRUE,
					   frequencyType = "rel")

  expect_error(oncoSimulIndiv(afe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = 20,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 5000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE),
	"If death is specified in the fitness effects, use Arb or Const model.")
})

test_that("testing output classes", {

  r <- data.frame(rfitness(2))

  r[, "Birth"] <- c("f_ - f_1 - f_2 - f_1_2",
                      "max(100*f_1, 10)",
                      "max(100*f_2, 10)",
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")

  r[, "Death"] <- c("f_ - f_1 - f_2 - f_1_2",
                      "max(100*f_1, 10)",
                      "max(100*f_2, 10)",
                      "max((200*(f_1 + f_2) + 50*f_1_2), 1)")


  afe <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  osi <- oncoSimulIndiv(afe,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 20,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 5000,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE)

  expect_identical(class(r), "data.frame")

  expect_identical(class(afe), c("fitnessEffects", "fitnessEffects_v3"))


  expect_identical(class(osi), c("oncosimul", "oncosimul2"))

  if(as.character(version$major) < 4) {
  expect_identical(class(osi$Genotypes), "matrix")
  } else {
      expect_identical(class(osi$Genotypes), c("matrix", "array"))
  }
})

test_that("testing performance", {

  r <- data.frame(rfitness(2))

  # Testing cases when totalPop should be 0
  r[, "Birth"] <- c(1, 1, 1, 1)

  r[, "Death"] <- c(1, 10, 10, 10)

  afe1 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = FALSE,
						   frequencyDependentDeath = FALSE,
						   deathSpec = TRUE)


  set.seed(1)

  null <- capture.output(osi1 <- oncoSimulIndiv(afe1,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = c(500, 500),
						initMutant = c("WT", "A"),
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  r[, "Birth"] <- c(1, 1, 1, 1)

  r[, "Death"] <- c("10*f_",
                    "10*f_1",
                    "50*f_2",
                    "200*(f_1 + f_2) + 50*f_1_2")



  afe2 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = FALSE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi2 <- oncoSimulIndiv(afe2,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))


  r[, "Birth"] <- c("f_", "f_1", "f_2", "f_1_2")

  r[, "Death"] <- c("10*f_",
                    "10*f_1",
                    "50*f_2",
                    "200*(f_1 + f_2) + 50*f_1_2")



  afe3 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi3 <- oncoSimulIndiv(afe3,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  r[, "Birth"] <- c("f_", "f_1", "f_2", "f_1_2")

  r[, "Death"] <- c(10, 10, 10, 10)



  afe4 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = FALSE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi4 <- oncoSimulIndiv(afe4,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  expect_output(print(osi1),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi1$NumClones, 2)

  expect_true(osi1$TotalPopSize == 0)

  expect_true(osi1$geneNames[2] %in% LETTERS)

  expect_output(print(osi2),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi2$NumClones, 1)

  expect_true(osi2$TotalPopSize == 0)


  expect_true(osi2$geneNames[2] %in% LETTERS)

  expect_true(osi2$FinalTime < 1)

  expect_output(print(osi3),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi3$NumClones, 1)

  expect_true(osi3$TotalPopSize == 0)


  expect_true(osi3$geneNames[2] %in% LETTERS)

  expect_true(osi3$FinalTime < 1)

  expect_output(print(osi4),
                "Individual OncoSimul trajectory", fixed = TRUE)


  expect_identical(osi4$NumClones, 1)

  expect_true(osi4$TotalPopSize == 0)


  expect_true(osi4$geneNames[2] %in% LETTERS)

  expect_true(osi4$FinalTime < 1)

  # Testing cases when the final totalPop is very big

  r[, "Death"] <- c(1, 1, 1, 1)

  r[, "Birth"] <- c(1, 15, 15, 15)

  afe5 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = FALSE,
						   frequencyDependentDeath = FALSE,
						   deathSpec = TRUE)


  set.seed(1)

  null <- capture.output(osi5 <- oncoSimulIndiv(afe5,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = c(500, 500),
						initMutant = c("WT", "A"),
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  r[, "Death"] <- c(1, 1, 1, 1)

  r[, "Birth"] <- c("10*f_",
                    "10*f_1",
                    "50*f_2",
                    "200*(f_1 + f_2) + 50*f_1_2")



  afe6 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = FALSE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi6 <- oncoSimulIndiv(afe6,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))


  r[, "Death"] <- c("f_", "f_1", "f_2", "f_1_2")

  r[, "Birth"] <- c("10*f_",
                    "10*f_1",
                    "50*f_2",
                    "200*(f_1 + f_2) + 50*f_1_2")



  afe7 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = TRUE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi7 <- oncoSimulIndiv(afe7,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  r[, "Death"] <- c("f_", "f_1", "f_2", "f_1_2")

  r[, "Birth"] <- c(10, 10, 10, 10)



  afe8 <- allFitnessEffects(genotFitness = r,
                           frequencyDependentBirth = FALSE,
						   frequencyDependentDeath = TRUE,
						   deathSpec = TRUE,
                           frequencyType = "rel")

  set.seed(1)

  null <- capture.output(osi8 <- oncoSimulIndiv(afe8,
                        model = "Arb",
                        onlyCancer = FALSE,
                        finalTime = 5000,
                        verbosity = 0,
                        mu = 1e-6,
                        initSize = 500,
                        keepPhylog = FALSE,
                        seed = NULL,
                        errorHitMaxTries = TRUE,
                        errorHitWallTime = TRUE))

  expect_output(print(osi5),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi5$NumClones, 3)

  expect_true(osi5$TotalPopSize > 1e8)


  expect_true(osi5$geneNames[2] %in% LETTERS)

  expect_output(print(osi6),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi6$NumClones, 3)

  expect_true(osi6$TotalPopSize > 1e8)


  expect_true(osi6$geneNames[2] %in% LETTERS)

  expect_true(osi6$FinalTime > 1)

  expect_output(print(osi7),
                "Individual OncoSimul trajectory", fixed = TRUE)

  ## Skip on kjohnson3, arm64
  if (Sys.getenv("R_PLATFORM") != "aarch64-apple-darwin20") {
    expect_identical(osi7$NumClones, 3)

    expect_true(osi7$TotalPopSize > 1e8)


    expect_true(osi7$geneNames[2] %in% LETTERS)

    expect_true(osi7$FinalTime > 1)
  }

  expect_output(print(osi8),
                "Individual OncoSimul trajectory", fixed = TRUE)

  expect_identical(osi8$NumClones, 3)

  expect_true(osi8$TotalPopSize > 1e8)


  expect_true(osi8$geneNames[2] %in% LETTERS)

  expect_true(osi8$FinalTime > 1)

})

set.seed(NULL)

cat(paste("\n Ending test.Z-oncoSimulIndivDeath at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

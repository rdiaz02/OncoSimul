# Non FDB example

genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c(1, 2, 3, 4),
                      stringsAsFactors = FALSE)

afe0 <- allFitnessEffects(genotFitness = genofit)

evalAllGenotypes(afe0, spPopSizes = c(5000, 2500, 2500, 500))

# Non FDB non FDD example
genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c(1, 2, 3, 4),
                      Death = c(1, 2, 3, 4),
                      stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genofit, deathSpec=TRUE)

evalAllGenotypes(afe, spPopSizes = c(5000, 2500, 2500, 500))

# Checking if the FDF (now FDB) works as it worked before
# 0's and 1's genotype

genofit <- data.frame(A = c(0, 1, 0, 1),
                     B = c(0, 0, 1, 1),
                     Birth = c("max(3, 2*f_)",
                                 "max(1.5, 3*(f_ + f_A))",
                                 "max(2, 3*(f_ + f_B))",
                                 "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                     stringsAsFactors = FALSE)

afe1 <- allFitnessEffects(genotFitness = genofit,
                  frequencyDependentBirth = TRUE,
                  frequencyType = "rel")

evalAllGenotypes(afe1, spPopSizes = c(5000, 2500, 2500, 500))
# print(afe)

# Letters genotype

genofit <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                      Birth = c("max(3, 2*f_)",
                                "max(1.5, 3*(f_ + f_A))",
                                "max(2, 3*(f_ + f_B))",
                                "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                      stringsAsFactors = FALSE)

afe2 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentBirth = TRUE,
                          frequencyType = "rel")

evalAllGenotypes(afe2, spPopSizes = c(5000, 2500, 2500, 500))

# print(afe2)


# Should trigger and error, because the genofit indicates 'rel' frequency and it
# is indicated 'abs'
afe3 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentBirth = TRUE,
                          frequencyType = "abs")



# Tests with FDB, but not FDD

genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c("max(3, 2*f_)",
                                "max(1.5, 3*(f_ + f_A))",
                                "max(2, 3*(f_ + f_B))",
                                "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                      Death = c(1, 2, 3, 4),
                      stringsAsFactors = FALSE)

afe3 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentBirth = TRUE,
                          frequencyType = "rel",
                          deathSpec=TRUE)

evalAllGenotypes(afe3, spPopSizes = c(5000, 2500, 2500, 500))

# Test with FDD, but not FDB

genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c(1, 2, 3, 4),
                      Death = c("max(3, 2*f_)",
                                "max(1.5, 3*(f_ + f_A))",
                                "max(2, 3*(f_ + f_B))",
                                "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                      stringsAsFactors = FALSE)

afe4 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentDeath = TRUE,
                          deathSpec=TRUE)

evalAllGenotypes(afe4, spPopSizes = c(5000, 2500, 2500, 500))
# Tests with both FDB and FDD


genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c("max(3, 2*f_)",
                                "max(1.5, 3*(f_ + f_A))",
                                "max(2, 3*(f_ + f_B))",
                                "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                      Death = c("1 + 1.5*f_",
                                "5 + 3*(f_A + f_B + f_A_B)",
                                "5 + 3*(f_A + f_B + f_A_B)",
                                "7 + 5*(f_A + f_B + f_A_B)"),
                      stringsAsFactors = FALSE)

afe5 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentBirth = TRUE,
                          frequencyDependentDeath = TRUE)

evalAllGenotypes(afe5, spPopSizes = c(5000, 2500, 2500, 500))

# Diferent frequencyType
genofit <- data.frame(A = c(0, 1, 0, 1),
                      B = c(0, 0, 1, 1),
                      Birth = c("max(3, 2*f_)",
                                "max(1.5, 3*(f_ + f_A))",
                                "max(2, 3*(f_ + f_B))",
                                "max(2, 5*f_ - 0.5*( f_A + f_B) + 15*f_A_B)"),
                      Death = c("1 + 1.5*n_",
                                "5 + 3*(n_A + n_B + n_A_B)",
                                "5 + 3*(n_A + n_B + n_A_B)",
                                "7 + 5*(n_A + n_B + n_A_B)"),
                      stringsAsFactors = FALSE)

afe6 <- allFitnessEffects(genotFitness = genofit,
                          frequencyDependentBirth = TRUE,
                          frequencyDependentDeath = TRUE)


osi <- oncoSimulIndiv(afe3,
                      model = "Arb",
                      onlyCancer = FALSE,
                      finalTime = 1,
                      mu = c("A" = 1e-6, B = 1e-8),
                      initMutant = c("WT", "A"),
                      initSize = c(5000, 1000),
                      keepPhylog = FALSE,
                      seed = NULL,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
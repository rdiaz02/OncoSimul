options(stringsAsFactors = FALSE)

r <- data.frame(Genotype = c("WT", "cucu", "cococ", "cucu, cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_A + f_B + f_A_B)",
                             "5 + 3*(f_A + f_B + f_A_B)",
                             "7 + 5*(f_A + f_B + f_A_B)"))

afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500, 500))

evalAllGenotypes(afe)

## pass spPopSizes to evalAllGenotyes, not in allFitnessEffects



r <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_A + f_B)",
                             "6 + 3*(f_A + f_B)"))
afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500))
evalAllGenotypes(afe)


r2 <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_1 + f_2)",
                             "6 + 3*(f_1 + f_2)"))
afe2 <- allFitnessEffects(genotFitness = r2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500))
evalAllGenotypes(afe2)





r2 <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1",
                             "1.5 + 1*(f_2)",
                             "2 + 1*(f_1)"))

afe3 <- allFitnessEffects(genotFitness = r2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

mt <- allMutatorEffects(epistasis = c("cucu" = 2,
                                      "cococ" = 4))


set.seed(1)
s1 <- oncoSimulIndiv(afe3, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-5,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s1, show = "genotypes")


set.seed(1)
s2 <- oncoSimulIndiv(afe3,
                     muEF = mt,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-5,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s2, show = "genotypes")


################


r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                 Fitness = c("1",
                             "1.4 + 1*(f_2)",
                             "1.4 + 1*(f_1)",
                             "1.6 + f_1 + f_2"))
afe4 <- allFitnessEffects(genotFitness = r4, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


set.seed(1)
s1 <- oncoSimulIndiv(afe4, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s1, show = "genotypes")


mt <- allMutatorEffects(epistasis = c("A" = 0.1,
                                      "B" = 10))
set.seed(1)
s2 <- oncoSimulIndiv(afe4,
                     muEF = mt,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s2, show = "genotypes")


plotClonePhylog(s1, keepEvents = TRUE)
plotClonePhylog(s2, keepEvents = TRUE)

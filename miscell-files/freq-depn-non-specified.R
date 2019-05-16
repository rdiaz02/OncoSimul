## Example that shows that non-specified are 0

options(stringsAsFactors = FALSE)

r1 <- data.frame(Genotype = c("WT", "A", "B", "C", "A, C", "A, B, C"), 
                 Fitness = c("1 + 1.5*f_",
                             "2 + 3 * f_B",
                             "3 + 2 * f_A",
                             "0",
                             "7",
                             "9"))
afe1 <- allFitnessEffects(genotFitness = r1, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


r2 <- data.frame(Genotype = c("WT", "A", "B", "A, C", "A, B, C"), 
                 Fitness = c("1 + 1.5*f_",
                             "2 + 3 * f_B",
                             "3 + 2 * f_A",
                             "7",
                             "9"))
afe2 <- allFitnessEffects(genotFitness = r2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


set.seed(1)
s1 <- oncoSimulIndiv(afe1, 
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


set.seed(1)
s2 <- oncoSimulIndiv(afe2, 
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



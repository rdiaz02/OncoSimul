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
                             "5 + 3*(f_A + f_B)"))

afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500))

evalAllGenotypes(afe)
